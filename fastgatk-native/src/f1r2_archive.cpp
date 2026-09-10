#include "fastgatk/somatic/f1r2_archive.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

namespace fastgatk::somatic {
namespace {

std::vector<std::string> all_contexts() {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    std::vector<std::string> result;
    result.reserve(64);
    for (const auto left : bases)
        for (const auto middle : bases)
            for (const auto right : bases)
                result.emplace_back(std::string{left, middle, right});
    return result;
}

char base_char(const int index) {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    return index >= 0 && index < 4 ? bases[static_cast<std::size_t>(index)] : 'N';
}

std::string url_encode(const std::string& value) {
    std::ostringstream out;
    out << std::uppercase << std::hex;
    for (const unsigned char c : value) {
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
            (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.' || c == '~')
            out << static_cast<char>(c);
        else
            out << '%' << std::setw(2) << std::setfill('0') << static_cast<int>(c)
                << std::setfill(' ');
    }
    return out.str();
}

std::uint64_t histogram_value(const std::vector<std::uint64_t>& histogram,
                              const int depth) {
    return depth >= 0 && static_cast<std::size_t>(depth) < histogram.size()
        ? histogram[static_cast<std::size_t>(depth)] : 0;
}

std::string histogram_content(const F1R2ArchiveSample& sample,
                              const std::vector<std::string>& contexts,
                              const bool reference_histogram) {
    std::ostringstream out;
    out << "## htsjdk.samtools.metrics.StringHeader\n# " << sample.sample
        << "\n\n\n## HISTOGRAM\tjava.lang.Integer\n";
    out << "depth";
    if (reference_histogram) {
        for (const auto& context : contexts) out << '\t' << context;
    } else {
        for (const auto& context : contexts) {
            const int middle = context[1] == 'A' ? 0 : context[1] == 'C' ? 1 :
                               context[1] == 'G' ? 2 : 3;
            for (int alt = 0; alt < 4; ++alt) {
                if (alt == middle) continue;
                out << '\t' << context << '_' << base_char(alt) << "_F1R2"
                    << '\t' << context << '_' << base_char(alt) << "_F2R1";
            }
        }
    }
    out << '\n';
    const int max_depth = std::max(1, sample.max_depth);
    for (int depth = 1; depth <= max_depth; ++depth) {
        out << depth;
        if (reference_histogram) {
            for (const auto& histogram : sample.ref_hist)
                out << '\t' << histogram_value(histogram, depth);
        } else {
            for (const auto& histogram : sample.alt_hist) {
                out << '\t' << histogram_value(histogram, depth);
                // alt_hist stores F1R2/F2R1 as adjacent columns.
                // The array is already laid out in that order.
            }
        }
        out << '\n';
    }
    return out.str();
}

std::string alt_table_content(const F1R2ArchiveSample& sample) {
    std::ostringstream out;
    out << "#<METADATA>SAMPLE=" << sample.sample << '\n'
        << "context\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt\n";
    for (const auto& row : sample.alt_rows) {
        out << row.context << '\t' << row.ref_count << '\t' << row.alt_count << '\t'
            << row.ref_f1r2 << '\t' << row.alt_f1r2 << '\t' << row.depth << '\t'
            << row.alt << '\n';
    }
    return out.str();
}

void append_tar_field(std::array<unsigned char, 512>& header,
                      const std::size_t offset, const std::size_t length,
                      const std::string& value) {
    const auto count = std::min(length, value.size());
    std::copy_n(value.data(), count, header.data() + offset);
}

void append_tar_octal(std::array<unsigned char, 512>& header,
                      const std::size_t offset, const std::size_t length,
                      const std::uint64_t value) {
    std::ostringstream octal;
    octal << std::oct << value;
    const auto text = octal.str();
    std::fill(header.begin() + offset, header.begin() + offset + length, '0');
    const auto begin = text.size() < length - 1 ? 0 : text.size() - (length - 1);
    std::copy(text.begin() + begin, text.end(),
              header.begin() + offset + (length - 1) - (text.size() - begin));
    header[offset + length - 1] = '\0';
}

void write_member(gzFile stream, const std::string& member,
                  const std::string& content) {
    std::array<unsigned char, 512> header{};
    append_tar_field(header, 0, 100, member);
    append_tar_field(header, 100, 8, "0000644\0");
    append_tar_field(header, 257, 6, "ustar\0");
    append_tar_field(header, 263, 2, "00");
    append_tar_octal(header, 124, 12, content.size());
    std::fill(header.begin() + 148, header.begin() + 156, ' ');
    header[156] = '0';
    std::uint32_t checksum = 0;
    for (const auto byte : header) checksum += byte;
    append_tar_octal(header, 148, 8, checksum);
    auto write = [&](const void* pointer, std::size_t bytes) {
        const auto* data = static_cast<const unsigned char*>(pointer);
        while (bytes != 0) {
            const auto chunk = static_cast<unsigned int>(std::min<std::size_t>(
                bytes, std::numeric_limits<unsigned int>::max()));
            if (gzwrite(stream, data, chunk) != static_cast<int>(chunk))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write F1R2 archive");
            data += chunk;
            bytes -= chunk;
        }
    };
    write(header.data(), header.size());
    write(content.data(), content.size());
    const std::array<unsigned char, 512> zeros{};
    const auto remainder = content.size() % 512;
    if (remainder != 0) write(zeros.data(), 512 - remainder);
}

}  // namespace

void write_f1r2_archive(const std::string& output,
                        const std::vector<F1R2ArchiveSample>& samples) {
    if (samples.empty())
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: no F1R2 samples");
    const auto contexts = all_contexts();
    gzFile stream = gzopen(output.c_str(), "wb");
    if (stream == nullptr)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create F1R2 archive: " + output);
    try {
        for (const auto& sample : samples) {
            const auto encoded = url_encode(sample.sample);
            write_member(stream, encoded + ".ref_histogram",
                         histogram_content(sample, contexts, true));
            write_member(stream, encoded + ".alt_histogram",
                         histogram_content(sample, contexts, false));
            write_member(stream, encoded + ".alt_table", alt_table_content(sample));
        }
        const std::array<unsigned char, 1024> zeros{};
        if (gzwrite(stream, zeros.data(), static_cast<unsigned int>(zeros.size())) !=
            static_cast<int>(zeros.size()))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize F1R2 archive");
        if (gzclose(stream) != Z_OK)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize F1R2 archive");
    } catch (...) {
        gzclose(stream);
        throw;
    }
}

}  // namespace fastgatk::somatic
