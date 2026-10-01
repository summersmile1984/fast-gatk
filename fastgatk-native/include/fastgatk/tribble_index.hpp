#pragma once

// Shared Tribble v3 index serialization.
//
// Extracted verbatim from index_feature_file_tool.cpp so that both
// IndexFeatureFile and CompareReferences write identical bytes.  GATK's
// CompareReferences creates its FULL_ALIGNMENT VCF through
// VariantContextWriterBuilder with an on-the-fly index, i.e. htsjdk's
// VCFIndexCreator + DynamicIndexCreator, which is the same serializer used for
// uncompressed VCF elsewhere -- with one addition: the VCF writer publishes a
// DICT:<contig>=<length> property for every contig in the VCF header before the
// four DynamicIndexCreator summary properties.

#include <sys/stat.h>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fastgatk::tribble {

inline constexpr int kTribbleMagic = 1480870228;
inline constexpr int kLinearIndexType = 1;
inline constexpr int kGvcfBinWidth = 128000;

struct LinearBlock {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

struct LinearContig {
    std::string name;
    std::vector<LinearBlock> blocks;
    int bin_width = kGvcfBinWidth;
    int longest_feature = 0;
    int n_features = 0;
};

struct LinearIndexStats {
    std::size_t records = 0;
    std::size_t contigs = 0;
    std::size_t blocks = 0;
    // DynamicIndexCreator publishes these values as v3 properties.  Keep the
    // summary separate from the linear layout so ordinary VCF can select the
    // same linear/interval-tree candidate as HTSJDK without carrying Java's
    // feature object graph into the native path.
    double feature_length_mean = 0.0;
    double feature_length_stddev = 0.0;
    double feature_length_variance = 0.0;
};

struct IntervalBlock {
    int start = 0;
    int end = 0;
    std::uint64_t file_start = 0;
    std::uint64_t file_end = 0;
};

struct IntervalContig {
    std::string name;
    std::vector<IntervalBlock> intervals;
};

struct FeatureSummary {
    std::size_t records = 0;
    std::uint64_t bases_seen = 0;
    int longest_feature = 0;
    // HTSJDK's RunningStat uses the usual online mean/M2 update.  The
    // DynamicIndexCreator pushes the running *maximum* feature length, not
    // the raw length, so retain that subtle behavior for property parity.
    double mean = 0.0;
    double m2 = 0.0;
    std::size_t stat_count = 0;
    std::vector<IntervalContig> interval_contigs;
};

inline void push_feature_stat(FeatureSummary& summary, double value) {
    ++summary.stat_count;
    const double delta = value - summary.mean;
    summary.mean += delta / static_cast<double>(summary.stat_count);
    const double delta2 = value - summary.mean;
    summary.m2 += delta * delta2;
}

inline double summary_variance(const FeatureSummary& summary) {
    // htsjdk.tribble.util.MathUtils.RunningStat reports the unbiased sample
    // variance (M2/(n-1)), not the population variance.  This is observable
    // in the BED fixture's FEATURE_LENGTH_STD_DEV and must be retained for
    // DynamicIndexCreator byte compatibility.
    return summary.stat_count <= 1 ? 0.0 :
        summary.m2 / static_cast<double>(summary.stat_count - 1);
}

inline double summary_stddev(const FeatureSummary& summary) {
    return std::sqrt(std::max(0.0, summary_variance(summary)));
}

inline std::string java_double_string(double value);
inline void optimize_linear_contig(LinearContig& contig) {
    // LinearIndex.optimize() adaptively merges adjacent bins until the
    // densest bin would contain more than MAX_FEATURES_PER_BIN (100) average
    // features, or the occupied-contig safety limit is reached.  This is
    // important for GVCFs: HTSJDK's advertised 128 kb starting bin is only a
    // starting point and the final .idx must carry the optimized width.
    constexpr double max_features_per_bin = 100.0;
    constexpr int max_bin_width = 1000000000;
    constexpr int max_occupied_bin_width = 1024000;
    // HTSJDK returns the last layout whose score still passes the threshold
    // (it saves "lastGood" *before* each merge and returns it when the merged
    // layout's score exceeds the threshold).  Mirror that so a merge that
    // would over-densify a bin is never emitted.
    LinearContig last_good = contig;
    while (contig.blocks.size() > 1) {
        std::uint64_t total_size = 0;
        for (const auto& block : contig.blocks) {
            const std::uint64_t size = block.end - block.start;
            total_size += size;
        }
        const double average_feature_size = static_cast<double>(total_size) /
            static_cast<double>(std::max(1, contig.n_features));
        double densest = -1.0;
        for (const auto& block : contig.blocks)
            densest = std::max(densest, static_cast<double>(block.end - block.start) /
                               average_feature_size);
        // HTSJDK's occupied-contig guard treats the configured maximum as a
        // terminal width (the next doubling is never attempted), so retain
        // the 1,024,000 bp layout observed in its serialized indices.  The
        // comparison is strict `>` to match htsjdk badBinWidth().
        const bool bad_width = contig.bin_width > max_bin_width || contig.bin_width < 0 ||
            (contig.n_features > 1 && contig.bin_width > max_occupied_bin_width);
        if (densest > max_features_per_bin || bad_width) {
            contig = last_good;
            return;
        }
        last_good = contig;
        std::vector<LinearBlock> merged;
        merged.reserve((contig.blocks.size() + 1) / 2);
        for (std::size_t index = 0; index < contig.blocks.size(); index += 2) {
            LinearBlock block = contig.blocks[index];
            if (index + 1 < contig.blocks.size())
                block.end += contig.blocks[index + 1].end - contig.blocks[index + 1].start;
            merged.push_back(block);
        }
        // HTSJDK's loop condition (`nBlocks != 1`) forbids collapsing an
        // occupied chromosome to a single block; restore the pre-merge layout.
        if (merged.size() == 1) {
            contig = last_good;
            return;
        }
        contig.blocks.swap(merged);
        if (contig.bin_width > max_bin_width / 2) contig.bin_width = max_bin_width + 1;
        else contig.bin_width *= 2;
    }
}

inline void write_i32(std::ofstream& output, std::int32_t value) {
    const std::uint32_t bits = static_cast<std::uint32_t>(value);
    for (int shift = 0; shift < 32; shift += 8)
        output.put(static_cast<char>((bits >> shift) & 0xffU));
}

inline void write_i64(std::ofstream& output, std::int64_t value) {
    const std::uint64_t bits = static_cast<std::uint64_t>(value);
    for (int shift = 0; shift < 64; shift += 8)
        output.put(static_cast<char>((bits >> shift) & 0xffU));
}

inline void write_tribble_string(std::ofstream& output, const std::string& value) {
    if (value.find('\0') != std::string::npos)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Tribble index strings cannot contain NUL");
    output.write(value.data(), static_cast<std::streamsize>(value.size()));
    output.put('\0');
}

inline std::string file_uri(const std::string& path) {
    const std::string absolute = std::filesystem::absolute(path).generic_string();
    std::ostringstream uri;
    uri << "file://";
    constexpr char hex[] = "0123456789ABCDEF";
    for (const unsigned char character : absolute) {
        const bool unreserved = (character >= 'a' && character <= 'z') ||
            (character >= 'A' && character <= 'Z') ||
            (character >= '0' && character <= '9') || character == '-' ||
            character == '_' || character == '.' || character == '~' || character == '/';
        if (unreserved) uri << static_cast<char>(character);
        else uri << '%' << hex[(character >> 4) & 0xf] << hex[character & 0xf];
    }
    return uri.str();
}

inline std::int64_t file_mtime_millis(const std::string& path) {
    struct stat status {};
    if (::stat(path.c_str(), &status) != 0)
        throw std::runtime_error("BAD_INPUT: cannot stat feature file: " + path);
    return static_cast<std::int64_t>(status.st_mtim.tv_sec) * 1000 +
        static_cast<std::int64_t>(status.st_mtim.tv_nsec / 1000000);
}
using TribbleProperties = std::vector<std::pair<std::string, std::string>>;

inline void write_tribble_header(std::ofstream& output, int index_type, const std::string& input_path,
                          std::uint64_t input_size, const TribbleProperties& properties);

// DynamicIndexCreator.finalizeIndex publishes exactly these four v3 properties
// on whichever candidate index it selects.
TribbleProperties dynamic_index_properties(const FeatureSummary& summary) {
    return TribbleProperties{
        {"FEATURE_LENGTH_MEAN", java_double_string(summary.mean)},
        {"FEATURE_LENGTH_STD_DEV", java_double_string(summary_stddev(summary))},
        {"MEAN_FEATURE_VARIANCE", java_double_string(summary_variance(summary))},
        {"FEATURE_COUNT", std::to_string(summary.records)}};
}
inline std::string java_double_string(double value) {
    // Reproduce java.lang.Double.toString exactly, because these strings are
    // serialized verbatim as Tribble v3 index properties (FEATURE_LENGTH_MEAN
    // and friends) and GATK compares them byte-for-byte.
    //
    // Java's contract is shortest-round-trip decimal digits, rendered as plain
    // decimal notation when 1e-3 <= |value| < 1e7 and as computerized
    // scientific notation otherwise.  Scientific exponents carry a '-' sign
    // but never a '+', and an integral mantissa always keeps its ".0"
    // ("1.0E7").  std::to_chars(shortest) supplies the same digits; only the
    // rendering rule and the exponent spelling have to be reconstructed.
    if (std::isnan(value)) return "NaN";
    if (std::isinf(value)) return value < 0 ? "-Infinity" : "Infinity";
    const bool negative = std::signbit(value);
    if (value == 0.0) return negative ? "-0.0" : "0.0";

    std::array<char, 64> buffer{};
    const auto converted = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                         std::fabs(value), std::chars_format::scientific);
    if (converted.ec != std::errc{})
        throw std::runtime_error("INTERNAL_ERROR: cannot format index property double");
    const std::string scientific(buffer.data(), converted.ptr);

    const auto marker = scientific.find('e');
    if (marker == std::string::npos)
        throw std::runtime_error("INTERNAL_ERROR: unexpected scientific double format");
    std::string digits;
    for (std::size_t index = 0; index < marker; ++index)
        if (scientific[index] != '.') digits.push_back(scientific[index]);
    const int exponent = std::stoi(scientific.substr(marker + 1));

    std::string text;
    if (exponent >= -3 && exponent < 7) {
        if (exponent >= 0) {
            const std::size_t whole = static_cast<std::size_t>(exponent) + 1;
            if (digits.size() <= whole) {
                text = digits + std::string(whole - digits.size(), '0') + ".0";
            } else {
                text = digits.substr(0, whole) + "." + digits.substr(whole);
            }
        } else {
            text = "0." + std::string(static_cast<std::size_t>(-exponent) - 1, '0') + digits;
        }
    } else {
        text = digits.substr(0, 1) + "." +
            (digits.size() > 1 ? digits.substr(1) : std::string("0")) +
            "E" + std::to_string(exponent);
    }
    return negative ? "-" + text : text;
}

inline void write_tribble_header(std::ofstream& output, int index_type, const std::string& input_path,
                          std::uint64_t input_size, const TribbleProperties& properties) {
    write_i32(output, kTribbleMagic);
    write_i32(output, index_type);
    write_i32(output, 3);  // AbstractIndex current version.
    write_tribble_string(output, file_uri(input_path));
    write_i64(output, static_cast<std::int64_t>(input_size));
    write_i64(output, file_mtime_millis(input_path));
    write_tribble_string(output, "");  // Java leaves MD5 unset.
    write_i32(output, 0);  // flags
    write_i32(output, static_cast<std::int32_t>(properties.size()));
    for (const auto& property : properties) {
        write_tribble_string(output, property.first);
        write_tribble_string(output, property.second);
    }
}

// IntervalTreeIndex is likewise codec-agnostic, so the same serializer serves
// VCF, BED and interval_list (index magic type 2).
inline LinearIndexStats write_interval_tree_index(const std::string& input_path,
                                           const std::string& index_path,
                                           std::uint64_t input_size,
                                           const FeatureSummary& summary,
                                           const TribbleProperties& properties =
                                               TribbleProperties{}) {
    std::ofstream output(index_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write Tribble index: " + index_path);
    const TribbleProperties effective = properties.empty()
        ? dynamic_index_properties(summary) : properties;
    write_tribble_header(output, 2, input_path, input_size, effective);
    write_i32(output, static_cast<std::int32_t>(summary.interval_contigs.size()));
    std::size_t interval_count = 0;

    // HTSJDK's IntervalTreeIndex does not serialize intervals in coordinate
    // order.  IntervalTree.getIntervals() is a pre-order walk of the
    // red-black tree built by IntervalTree.insert(): equal starts are sent to
    // the left child, then the usual CLRS insert-fixup rotations are applied.
    // Reproduce that tree shape here so native dense indexes can be compared
    // byte-for-byte with DynamicIndexCreator, rather than only being query
    // compatible.  The min/max augmentation is used by HTSJDK for queries but
    // is not serialized; the topology and insertion comparator determine the
    // emitted order.
    struct TreeNode {
        IntervalBlock interval;
        int parent = -1;
        int left = -1;
        int right = -1;
        bool red = true;
    };

    auto emit_node = [&](const TreeNode& node) {
        if (node.interval.end < node.interval.start ||
            node.interval.file_end < node.interval.file_start ||
            node.interval.file_end - node.interval.file_start >
                static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: invalid interval block while writing Tribble index");
        }
        write_i32(output, node.interval.start);
        write_i32(output, node.interval.end);
        write_i64(output, static_cast<std::int64_t>(node.interval.file_start));
        write_i32(output, static_cast<std::int32_t>(
            node.interval.file_end - node.interval.file_start));
        ++interval_count;
    };

    for (const auto& contig : summary.interval_contigs) {
        write_tribble_string(output, contig.name);
        write_i32(output, static_cast<std::int32_t>(contig.intervals.size()));
        std::vector<TreeNode> nodes;
        nodes.reserve(contig.intervals.size());
        int root = -1;

        auto is_red = [&](const int index) {
            return index >= 0 && nodes[static_cast<std::size_t>(index)].red;
        };
        auto left_rotate = [&](const int x) {
            const int y = nodes[static_cast<std::size_t>(x)].right;
            if (y < 0) throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: malformed interval tree rotation");
            nodes[static_cast<std::size_t>(x)].right =
                nodes[static_cast<std::size_t>(y)].left;
            if (nodes[static_cast<std::size_t>(y)].left >= 0)
                nodes[static_cast<std::size_t>(nodes[static_cast<std::size_t>(y)].left)].parent = x;
            nodes[static_cast<std::size_t>(y)].parent =
                nodes[static_cast<std::size_t>(x)].parent;
            if (nodes[static_cast<std::size_t>(x)].parent < 0) {
                root = y;
            } else if (x == nodes[static_cast<std::size_t>(
                           nodes[static_cast<std::size_t>(x)].parent)].left) {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].left = y;
            } else {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].right = y;
            }
            nodes[static_cast<std::size_t>(y)].left = x;
            nodes[static_cast<std::size_t>(x)].parent = y;
        };
        auto right_rotate = [&](const int x) {
            const int y = nodes[static_cast<std::size_t>(x)].left;
            if (y < 0) throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: malformed interval tree rotation");
            nodes[static_cast<std::size_t>(x)].left =
                nodes[static_cast<std::size_t>(y)].right;
            if (nodes[static_cast<std::size_t>(y)].right >= 0)
                nodes[static_cast<std::size_t>(nodes[static_cast<std::size_t>(y)].right)].parent = x;
            nodes[static_cast<std::size_t>(y)].parent =
                nodes[static_cast<std::size_t>(x)].parent;
            if (nodes[static_cast<std::size_t>(x)].parent < 0) {
                root = y;
            } else if (x == nodes[static_cast<std::size_t>(
                           nodes[static_cast<std::size_t>(x)].parent)].right) {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].right = y;
            } else {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].left = y;
            }
            nodes[static_cast<std::size_t>(y)].right = x;
            nodes[static_cast<std::size_t>(x)].parent = y;
        };
        auto insert_fixup = [&](int z) {
            while (z != root && is_red(nodes[static_cast<std::size_t>(z)].parent)) {
                const int parent = nodes[static_cast<std::size_t>(z)].parent;
                const int grand = nodes[static_cast<std::size_t>(parent)].parent;
                if (parent == nodes[static_cast<std::size_t>(grand)].left) {
                    const int uncle = nodes[static_cast<std::size_t>(grand)].right;
                    if (is_red(uncle)) {
                        nodes[static_cast<std::size_t>(parent)].red = false;
                        nodes[static_cast<std::size_t>(uncle)].red = false;
                        nodes[static_cast<std::size_t>(grand)].red = true;
                        z = grand;
                    } else {
                        if (z == nodes[static_cast<std::size_t>(parent)].right) {
                            z = parent;
                            left_rotate(z);
                        }
                        const int fixed_parent = nodes[static_cast<std::size_t>(z)].parent;
                        const int fixed_grand = nodes[static_cast<std::size_t>(fixed_parent)].parent;
                        nodes[static_cast<std::size_t>(fixed_parent)].red = false;
                        nodes[static_cast<std::size_t>(fixed_grand)].red = true;
                        right_rotate(fixed_grand);
                    }
                } else {
                    const int uncle = nodes[static_cast<std::size_t>(grand)].left;
                    if (is_red(uncle)) {
                        nodes[static_cast<std::size_t>(parent)].red = false;
                        nodes[static_cast<std::size_t>(uncle)].red = false;
                        nodes[static_cast<std::size_t>(grand)].red = true;
                        z = grand;
                    } else {
                        if (z == nodes[static_cast<std::size_t>(parent)].left) {
                            z = parent;
                            right_rotate(z);
                        }
                        const int fixed_parent = nodes[static_cast<std::size_t>(z)].parent;
                        const int fixed_grand = nodes[static_cast<std::size_t>(fixed_parent)].parent;
                        nodes[static_cast<std::size_t>(fixed_parent)].red = false;
                        nodes[static_cast<std::size_t>(fixed_grand)].red = true;
                        left_rotate(fixed_grand);
                    }
                }
            }
            if (root >= 0) nodes[static_cast<std::size_t>(root)].red = false;
        };

        for (const auto& interval : contig.intervals) {
            int parent = -1;
            int cursor = root;
            while (cursor >= 0) {
                parent = cursor;
                // HTSJDK's treeInsert uses <= for the left branch, so equal
                // starts are deliberately not made stable on the right.
                cursor = interval.start <= nodes[static_cast<std::size_t>(cursor)].interval.start
                    ? nodes[static_cast<std::size_t>(cursor)].left
                    : nodes[static_cast<std::size_t>(cursor)].right;
            }
            const int inserted = static_cast<int>(nodes.size());
            nodes.push_back(TreeNode{interval, parent, -1, -1, true});
            if (parent < 0) root = inserted;
            else if (interval.start <= nodes[static_cast<std::size_t>(parent)].interval.start)
                nodes[static_cast<std::size_t>(parent)].left = inserted;
            else nodes[static_cast<std::size_t>(parent)].right = inserted;
            insert_fixup(inserted);
        }

        std::function<void(int)> preorder = [&](const int index) {
            if (index < 0) return;
            emit_node(nodes[static_cast<std::size_t>(index)]);
            preorder(nodes[static_cast<std::size_t>(index)].left);
            preorder(nodes[static_cast<std::size_t>(index)].right);
        };
        preorder(root);
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing Tribble interval index: " + index_path);
    return {summary.records, summary.interval_contigs.size(), interval_count,
            summary.mean, summary_stddev(summary), summary_variance(summary)};
}

// DynamicIndexCreator.finalizeIndex scores the linear candidate with
// linearScore = binSize * (featureCount / basesSeen) * ceil(longestFeature /
// binSize) and, under FOR_SEEK_TIME, picks the 2000 bp linear index only when
// that score stays below the interval tree's 75-features-per-interval budget.
// An empty file has no linear candidate at all.
inline bool prefers_linear_index(const FeatureSummary& summary) {
    if (summary.records == 0) return false;
    const double density = summary.bases_seen == 0 ? 0.0 :
        static_cast<double>(summary.records) / static_cast<double>(summary.bases_seen);
    const double linear_score = 2000.0 * density *
        std::ceil(static_cast<double>(summary.longest_feature) / 2000.0);
    return linear_score < 75.0;
}

// ---------------------------------------------------------------------------
// Non-streaming convenience layer.
//
// IndexFeatureFile streams its input so that huge GVCFs never materialise a
// feature list; CompareReferences already holds every record it emitted, so it
// uses these helpers to build the identical DynamicIndexCreator summary and
// LinearIndexCreator layout from an in-memory feature list.  The rules below
// mirror summarize_text_features()/write_linear_text_index() exactly.
// ---------------------------------------------------------------------------

struct ScannedFeature {
    std::string contig;
    int start = 0;
    int end = 0;
    std::uint64_t offset = 0;
};

struct TextIndexPlan {
    FeatureSummary summary;
    std::vector<LinearContig> linear_contigs;
};

inline TextIndexPlan plan_text_index(const std::vector<ScannedFeature>& features,
                                     std::uint64_t input_size, int initial_bin_width,
                                     bool want_linear) {
    TextIndexPlan plan;
    FeatureSummary& summary = plan.summary;
    std::string last_contig;
    std::string current_contig;
    int last_start = 0;
    std::size_t records_in_interval = 0;
    std::map<std::string, std::size_t> seen_contigs;
    int linear_last_start = 0;
    for (const auto& feature : features) {
        if (summary.records != 0) {
            if (feature.contig == last_contig && feature.start >= last_start)
                summary.bases_seen += static_cast<std::uint64_t>(feature.start - last_start);
            else
                summary.bases_seen += static_cast<std::uint64_t>(feature.start);
        } else {
            summary.bases_seen = static_cast<std::uint64_t>(feature.start);
        }
        if (feature.contig == last_contig && feature.start < last_start)
            throw std::runtime_error("BAD_INPUT: feature records are not sorted by coordinate");

        const int length = feature.end - feature.start + 1;
        summary.longest_feature = std::max(summary.longest_feature, length);
        push_feature_stat(summary, static_cast<double>(summary.longest_feature));
        ++summary.records;
        last_contig = feature.contig;
        last_start = feature.start;

        if (feature.contig != current_contig) {
            if (!current_contig.empty() && !summary.interval_contigs.empty() &&
                !summary.interval_contigs.back().intervals.empty())
                summary.interval_contigs.back().intervals.back().file_end = feature.offset;
            for (const auto& existing : summary.interval_contigs)
                if (existing.name == feature.contig)
                    throw std::runtime_error(
                        "BAD_INPUT: feature contigs are not grouped in coordinate order: " + feature.contig);
            summary.interval_contigs.push_back(IntervalContig{feature.contig, {}});
            current_contig = feature.contig;
            records_in_interval = 0;
        }
        auto& intervals = summary.interval_contigs.back().intervals;
        if (intervals.empty() || records_in_interval >= 75) {
            if (!intervals.empty()) intervals.back().file_end = feature.offset;
            intervals.push_back(IntervalBlock{feature.start, feature.end, feature.offset, 0});
            records_in_interval = 0;
        }
        intervals.back().end = std::max(intervals.back().end, feature.end);
        ++records_in_interval;

        if (!want_linear) continue;
        if (plan.linear_contigs.empty() || plan.linear_contigs.back().name != feature.contig) {
            if (seen_contigs.find(feature.contig) != seen_contigs.end())
                throw std::runtime_error(
                    "BAD_INPUT: feature contigs are not grouped in coordinate order: " + feature.contig);
            seen_contigs.emplace(feature.contig, plan.linear_contigs.size());
            plan.linear_contigs.push_back(
                LinearContig{feature.contig, {{feature.offset, 0}}, initial_bin_width, 0, 0});
            linear_last_start = 0;
        }
        if (feature.start < linear_last_start)
            throw std::runtime_error("BAD_INPUT: feature records are not sorted by coordinate");
        linear_last_start = feature.start;
        auto& contig = plan.linear_contigs.back();
        while (static_cast<std::int64_t>(feature.start) >
               static_cast<std::int64_t>(contig.blocks.size()) * initial_bin_width)
            contig.blocks.push_back({feature.offset, 0});
        contig.longest_feature = std::max(contig.longest_feature, length);
        ++contig.n_features;
    }
    if (!summary.interval_contigs.empty() && !summary.interval_contigs.back().intervals.empty())
        summary.interval_contigs.back().intervals.back().file_end = input_size;
    for (std::size_t index = 0; index < plan.linear_contigs.size(); ++index) {
        auto& contig = plan.linear_contigs[index];
        const std::uint64_t contig_end = index + 1 < plan.linear_contigs.size()
            ? plan.linear_contigs[index + 1].blocks.front().start : input_size;
        for (std::size_t block = 0; block < contig.blocks.size(); ++block) {
            contig.blocks[block].end = block + 1 < contig.blocks.size()
                ? contig.blocks[block + 1].start : contig_end;
            if (contig.blocks[block].end < contig.blocks[block].start)
                throw std::runtime_error("BAD_INPUT: invalid feature block offsets");
        }
        optimize_linear_contig(contig);
    }
    return plan;
}

inline LinearIndexStats write_linear_index_body(const std::string& input_path,
                                                const std::string& index_path,
                                                std::uint64_t input_size, int index_type,
                                                const TribbleProperties& properties,
                                                const std::vector<LinearContig>& contigs) {
    std::ofstream output(index_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write Tribble index: " + index_path);
    write_tribble_header(output, index_type, input_path, input_size, properties);
    write_i32(output, static_cast<std::int32_t>(contigs.size()));
    std::size_t block_count = 0;
    std::size_t record_count = 0;
    for (const auto& contig : contigs) {
        write_tribble_string(output, contig.name);
        write_i32(output, contig.bin_width);
        write_i32(output, static_cast<std::int32_t>(contig.blocks.size()));
        write_i32(output, contig.longest_feature);
        write_i32(output, 0);  // OLD_V3_INDEX=false
        write_i32(output, contig.n_features);
        std::uint64_t previous_end = 0;
        for (const auto& block : contig.blocks) {
            write_i64(output, static_cast<std::int64_t>(block.start));
            previous_end = block.end;
        }
        write_i64(output, static_cast<std::int64_t>(previous_end));
        block_count += contig.blocks.size();
        record_count += static_cast<std::size_t>(contig.n_features);
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing Tribble index: " + index_path);
    return LinearIndexStats{record_count, contigs.size(), block_count};
}

}  // namespace fastgatk::tribble
