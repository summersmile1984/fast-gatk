#include "fastgatk/kernels/smith_waterman.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fastgatk::kernels {

namespace {

constexpr int kMatrixMinCutoff = static_cast<int>(-1.0e8);
constexpr int kLowInitValue = std::numeric_limits<int>::min() / 2;

enum class TraceState : std::uint8_t { Match, Insertion, Deletion, Clip };

struct Matrix {
    std::size_t rows = 0;     // reference length + 1
    std::size_t columns = 0; // alternate/read length + 1
    std::vector<int> scores;
    std::vector<int> backtrack;

    int& score(std::size_t row, std::size_t column) {
        return scores[row * columns + column];
    }
    int score(std::size_t row, std::size_t column) const {
        return scores[row * columns + column];
    }
    int& trace(std::size_t row, std::size_t column) {
        return backtrack[row * columns + column];
    }
    int trace(std::size_t row, std::size_t column) const {
        return backtrack[row * columns + column];
    }
};

bool has_gap_initialized_edges(SmithWatermanOverhangStrategy strategy) {
    return strategy == SmithWatermanOverhangStrategy::Indel ||
           strategy == SmithWatermanOverhangStrategy::LeadingIndel;
}

Matrix calculate_matrix(const std::uint8_t* reference, std::size_t reference_length,
                        const std::uint8_t* alternate, std::size_t alternate_length,
                        SmithWatermanParameters parameters,
                        SmithWatermanOverhangStrategy strategy) {
    Matrix matrix;
    matrix.rows = reference_length + 1;
    matrix.columns = alternate_length + 1;
    matrix.scores.assign(matrix.rows * matrix.columns, 0);
    matrix.backtrack.assign(matrix.rows * matrix.columns, 0);

    std::vector<int> best_gap_vertical(alternate_length + 1, kLowInitValue);
    std::vector<int> vertical_gap_size(alternate_length + 1, 0);
    std::vector<int> best_gap_horizontal(reference_length + 1, kLowInitValue);
    std::vector<int> horizontal_gap_size(reference_length + 1, 0);

    if (has_gap_initialized_edges(strategy)) {
        if (alternate_length >= 1) matrix.score(0, 1) = parameters.gap_open;
        for (std::size_t column = 2; column <= alternate_length; ++column)
            matrix.score(0, column) = matrix.score(0, column - 1) + parameters.gap_extend;
        if (reference_length >= 1) matrix.score(1, 0) = parameters.gap_open;
        for (std::size_t row = 2; row <= reference_length; ++row)
            matrix.score(row, 0) = matrix.score(row - 1, 0) + parameters.gap_extend;
    }

    for (std::size_t row = 1; row <= reference_length; ++row) {
        for (std::size_t column = 1; column <= alternate_length; ++column) {
            const int diagonal = matrix.score(row - 1, column - 1) +
                (reference[row - 1] == alternate[column - 1]
                     ? parameters.match : parameters.mismatch);

            const int previous_vertical = matrix.score(row - 1, column) + parameters.gap_open;
            best_gap_vertical[column] += parameters.gap_extend;
            if (previous_vertical > best_gap_vertical[column]) {
                best_gap_vertical[column] = previous_vertical;
                vertical_gap_size[column] = 1;
            } else {
                ++vertical_gap_size[column];
            }

            const int previous_horizontal = matrix.score(row, column - 1) + parameters.gap_open;
            best_gap_horizontal[row] += parameters.gap_extend;
            if (previous_horizontal > best_gap_horizontal[row]) {
                best_gap_horizontal[row] = previous_horizontal;
                horizontal_gap_size[row] = 1;
            } else {
                ++horizontal_gap_size[row];
            }

            const int down = best_gap_vertical[column];
            const int right = best_gap_horizontal[row];
            // GATK's Java aligner uses diagonal >= right/down, then right >=
            // down.  Keeping this order is important for homopolymer CIGARs.
            if (diagonal >= down && diagonal >= right) {
                matrix.score(row, column) = std::max(kMatrixMinCutoff, diagonal);
                matrix.trace(row, column) = 0;
            } else if (right >= down) {
                matrix.score(row, column) = std::max(kMatrixMinCutoff, right);
                matrix.trace(row, column) = -horizontal_gap_size[row];
            } else {
                matrix.score(row, column) = std::max(kMatrixMinCutoff, down);
                matrix.trace(row, column) = vertical_gap_size[column];
            }
        }
    }
    return matrix;
}

struct Endpoint {
    std::size_t reference = 0;
    std::size_t alternate = 0;
    int score = std::numeric_limits<int>::min();
};

Endpoint choose_endpoint(const Matrix& matrix, SmithWatermanOverhangStrategy strategy) {
    const auto reference_length = matrix.rows - 1;
    const auto alternate_length = matrix.columns - 1;
    if (strategy == SmithWatermanOverhangStrategy::Indel)
        return {reference_length, alternate_length,
                matrix.score(reference_length, alternate_length)};

    Endpoint endpoint{0, alternate_length, std::numeric_limits<int>::min()};
    // First inspect the rightmost column.  >= intentionally selects the
    // endpoint closest to the diagonal when scores tie.
    for (std::size_t row = 1; row <= reference_length; ++row) {
        const int score = matrix.score(row, alternate_length);
        if (score >= endpoint.score)
            endpoint = {row, alternate_length, score};
    }
    if (strategy != SmithWatermanOverhangStrategy::LeadingIndel) {
        for (std::size_t column = 1; column <= alternate_length; ++column) {
            const int score = matrix.score(reference_length, column);
            const auto current_distance = endpoint.reference > endpoint.alternate
                ? endpoint.reference - endpoint.alternate
                : endpoint.alternate - endpoint.reference;
            const auto candidate_distance = reference_length > column
                ? reference_length - column : column - reference_length;
            if (score > endpoint.score ||
                (score == endpoint.score && candidate_distance < current_distance))
                endpoint = {reference_length, column, score};
        }
    }
    return endpoint;
}

std::string compact_cigar(const std::vector<std::pair<TraceState, std::size_t>>& segments) {
    std::string cigar;
    TraceState previous = TraceState::Clip;
    std::size_t run = 0;
    auto flush = [&]() {
        if (run == 0) return;
        cigar += std::to_string(run);
        switch (previous) {
            case TraceState::Match: cigar.push_back('M'); break;
            case TraceState::Insertion: cigar.push_back('I'); break;
            case TraceState::Deletion: cigar.push_back('D'); break;
            case TraceState::Clip: cigar.push_back('S'); break;
        }
        run = 0;
    };
    for (const auto& [state, length] : segments) {
        if (length == 0) continue;
        if (run != 0 && state != previous) flush();
        previous = state;
        run += length;
    }
    flush();
    return cigar;
}

std::vector<std::pair<TraceState, std::size_t>> traceback(
    const Matrix& matrix, Endpoint endpoint,
    SmithWatermanOverhangStrategy strategy, std::size_t& reference_start,
    std::size_t& alternate_start) {
    std::size_t row = endpoint.reference;
    std::size_t column = endpoint.alternate;
    std::vector<std::pair<TraceState, std::size_t>> reverse_segments;
    const auto trailing_overhang = (matrix.columns - 1) - endpoint.alternate;
    // Java's calculateCigar carries a bottom-row trailing overhang in
    // segment_length.  SOFTCLIP emits it as a trailing S; the other modes
    // fold it into the aligned segment according to their strategy.
    std::size_t segment_length = trailing_overhang;
    if (strategy == SmithWatermanOverhangStrategy::Softclip && trailing_overhang != 0) {
        reverse_segments.emplace_back(TraceState::Clip, trailing_overhang);
        segment_length = 0;
    }
    TraceState state = TraceState::Match;
    do {
        const int backtrack = matrix.trace(row, column);
        TraceState next_state = TraceState::Match;
        std::size_t step_length = 1;
        if (backtrack > 0) {
            next_state = TraceState::Deletion;
            step_length = static_cast<std::size_t>(backtrack);
        } else if (backtrack < 0) {
            next_state = TraceState::Insertion;
            step_length = static_cast<std::size_t>(-backtrack);
        }

        if (next_state == state) {
            segment_length += step_length;
        } else {
            if (segment_length != 0) reverse_segments.emplace_back(state, segment_length);
            segment_length = step_length;
            state = next_state;
        }
        switch (next_state) {
            case TraceState::Match: --row; --column; break;
            case TraceState::Insertion: column -= step_length; break;
            case TraceState::Deletion: row -= step_length; break;
            case TraceState::Clip: break;
        }
    } while (row > 0 && column > 0);

    if (strategy == SmithWatermanOverhangStrategy::Softclip) {
        reverse_segments.emplace_back(state, segment_length);
        if (column > 0) reverse_segments.emplace_back(TraceState::Clip, column);
        reference_start = row;
        alternate_start = column;
    } else if (strategy == SmithWatermanOverhangStrategy::Ignore) {
        reverse_segments.emplace_back(state, segment_length + column);
        reference_start = row;
        alternate_start = column;
    } else {
        reverse_segments.emplace_back(state, segment_length);
        if (row > 0) reverse_segments.emplace_back(TraceState::Deletion, row);
        else if (column > 0) reverse_segments.emplace_back(TraceState::Insertion, column);
        reference_start = 0;
        alternate_start = 0;
    }
    std::reverse(reverse_segments.begin(), reverse_segments.end());
    return reverse_segments;
}

std::size_t last_substring(const std::uint8_t* reference, std::size_t reference_length,
                           const std::uint8_t* alternate, std::size_t alternate_length) {
    if (alternate_length > reference_length) return std::string::npos;
    for (std::size_t start = reference_length - alternate_length + 1; start-- > 0;) {
        bool equal = true;
        for (std::size_t offset = 0; offset < alternate_length; ++offset) {
            if (reference[start + offset] != alternate[offset]) {
                equal = false;
                break;
            }
        }
        if (equal) return start;
        if (start == 0) break;
    }
    return std::string::npos;
}

}  // namespace

int smith_waterman_score_reference(const std::uint8_t* read, std::size_t read_length,
                                   const std::uint8_t* reference, std::size_t reference_length,
                                   SmithWatermanParameters parameters,
                                   SmithWatermanOverhangStrategy overhang) {
    if ((read_length != 0 && read == nullptr) ||
        (reference_length != 0 && reference == nullptr))
        throw std::invalid_argument("Smith-Waterman input pointer is null");
    if (read_length == 0 || reference_length == 0)
        throw std::invalid_argument("Smith-Waterman sequences must be non-empty");
    if ((overhang == SmithWatermanOverhangStrategy::Softclip ||
         overhang == SmithWatermanOverhangStrategy::Ignore) &&
        last_substring(reference, reference_length, read, read_length) != std::string::npos)
        return static_cast<int>(read_length) * parameters.match;
    const auto matrix = calculate_matrix(reference, reference_length, read, read_length,
                                         parameters, overhang);
    return choose_endpoint(matrix, overhang).score;
}

SmithWatermanAlignment smith_waterman_align_reference(
    const std::uint8_t* read, std::size_t read_length,
    const std::uint8_t* reference, std::size_t reference_length,
    SmithWatermanParameters parameters, SmithWatermanOverhangStrategy overhang) {
    if ((read_length != 0 && read == nullptr) ||
        (reference_length != 0 && reference == nullptr))
        throw std::invalid_argument("Smith-Waterman input pointer is null");
    if (read_length == 0 || reference_length == 0)
        throw std::invalid_argument("Smith-Waterman sequences must be non-empty");

    SmithWatermanAlignment result;
    const auto exact = (overhang == SmithWatermanOverhangStrategy::Softclip ||
                        overhang == SmithWatermanOverhangStrategy::Ignore)
        ? last_substring(reference, reference_length, read, read_length)
        : std::string::npos;
    if (exact != std::string::npos) {
        result.score = static_cast<int>(read_length) * parameters.match;
        result.read_start = 0;
        result.read_end = read_length;
        result.reference_start = exact;
        result.reference_end = exact + read_length;
        result.alignment_offset = static_cast<int>(exact);
        result.cigar = std::to_string(read_length) + "M";
        return result;
    }

    const auto matrix = calculate_matrix(reference, reference_length, read, read_length,
                                         parameters, overhang);
    const auto endpoint = choose_endpoint(matrix, overhang);
    std::size_t reference_start = 0, alternate_start = 0;
    const auto segments = traceback(matrix, endpoint, overhang, reference_start, alternate_start);
    result.score = endpoint.score;
    result.read_start = alternate_start;
    result.read_end = endpoint.alternate;
    result.reference_start = reference_start;
    result.reference_end = endpoint.reference;
    if (overhang == SmithWatermanOverhangStrategy::Softclip)
        result.alignment_offset = static_cast<int>(reference_start);
    else if (overhang == SmithWatermanOverhangStrategy::Ignore)
        result.alignment_offset = static_cast<int>(reference_start) -
            static_cast<int>(alternate_start);
    result.cigar = compact_cigar(segments);
    return result;
}

}  // namespace fastgatk::kernels
