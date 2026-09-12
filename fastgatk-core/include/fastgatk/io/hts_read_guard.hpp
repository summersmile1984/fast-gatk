#pragma once

// One place where "htslib could not read this record" stops being mistaken for
// "htslib reached the end of the input".
//
// Why this header exists
// ----------------------
// Every htslib record reader reports a *clean end of input* and a *failed
// read* through two different negative return values:
//
//   bcf_read()    "0 on success; -1 on end of file; < -1 on critical error"
//                 (htslib/vcf.h; bcf_read -> vcf_read -> vcf_parse, whose
//                 ``err:`` label returns -2 for a record it rejected)
//   sam_read1()   "-1 at end of file; < -1 on failure" (htslib/sam.h)
//   hts_getline() "-1 on end-of-file; <= -2 on error" (htslib/hts.h; hts.c
//                 sets -2 when the underlying read itself failed)
//
// htslib states the rule itself where it consumes its own reader: in
// ``hts_readlist`` (hts.c) the loop is ``while ((ret = bgzf_getline(...)) >= 0)``
// followed by ``if (ret < -1) // Read error`` -- the two negatives are never
// collapsed.
//
// A loop written as ``while (bcf_read(fp, hdr, rec) == 0)`` collapses them: the
// FIRST record htslib could not read ends the traversal exactly as if the input
// had ended, so that record *and every record after it* are discarded while the
// process still exits 0.  The measured consequence on a truncated ``.vcf.gz`` is
// exit 0 with 19322 of 20001 records written and nothing on stdout to say so
// (see .diag/round-bcfread-audit.md).
//
// These wrappers perform the read and raise BAD_INPUT on the error, so a loop
// condition can only ever observe 0 (a record was read) or -1 (clean end of
// input).  They are deliberately thin: the loop body, the error convention
// (``throw std::runtime_error("BAD_INPUT: ...")`` -> exit 2) and every tool's
// behaviour on well-formed input are untouched.

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <stdexcept>
#include <string>

namespace fastgatk::io {

// htslib has no accessor for a per-file error string, so the diagnostic carries
// the numeric status plus the path; htslib has already logged its own
// `[E::vcf_parse_*]` / `[E::bgzf_read]` line naming the offending record.
[[noreturn]] inline void throw_read_failure(const char* what, int status,
                                            const std::string& path) {
    throw std::runtime_error(std::string("BAD_INPUT: ") + what +
                             " (htslib status " + std::to_string(status) + "): " + path);
}

// bcf_read / vcf_read.  Returns 0 or -1 only; a rejected record is reported.
inline int read_variant_record(htsFile* file, const bcf_hdr_t* header,
                               bcf1_t* record, const std::string& path) {
    const int status = bcf_read(file, header, record);
    if (status < -1)
        throw_read_failure("malformed or truncated VCF/BCF record", status, path);
    return status;
}

// sam_read1.  Returns 0 or -1 only; a corrupt record is reported.
inline int read_alignment_record(samFile* file, sam_hdr_t* header,
                                 bam1_t* record, const std::string& path) {
    const int status = sam_read1(file, header, record);
    if (status < -1)
        throw_read_failure("malformed or truncated BAM/CRAM record", status, path);
    return status;
}

// sam_itr_next: -1 means this iterator is exhausted, < -1 means the read failed.
inline int read_indexed_alignment_record(samFile* file, hts_itr_t* iterator,
                                         bam1_t* record, const std::string& path) {
    const int status = sam_itr_next(file, iterator, record);
    if (status < -1)
        throw_read_failure("malformed or truncated indexed BAM/CRAM record", status, path);
    return status;
}

// hts_getline: -1 means end of file, <= -2 means the stream itself failed (a
// truncated gzip/deflate block, for instance).  hts_getline accepts '\n' or
// KS_SEP_LINE (htslib/kseq.h), which have the same value for this purpose.
inline int read_text_line(htsFile* file, kstring_t* line, const std::string& path) {
    const int status = hts_getline(file, '\n', line);
    if (status < -1)
        throw_read_failure("malformed or truncated text record", status, path);
    return status;
}

}  // namespace fastgatk::io
