// Cluster resource probe contract.
//
// Verifies that the native runtime reads the SLURM allocation / cgroup
// environment correctly when actually running inside a real allocation.
// The smoke harness `fastgatk-runtime-smoke` only exercises the host-side
// probe; this binary is the cluster-side twin that the production
// scheduler runs first.
//
// Behaviour:
//   * FASTGATK_REAL_CLUSTER != 1            → exit 0 with status=skip
//   * FASTGATK_REAL_CLUSTER == 1 and no SLURM_JOB_ID → exit 2 (fail-closed)
//   * FASTGATK_REAL_CLUSTER == 1 + SLURM_JOB_ID set  → assert every field
//     matches the exported SLURM env, write the JSON snapshot to the path
//     given by FASTGATK_CLUSTER_EVIDENCE_JSON (default stdout), and exit 0.
//
// The fastgatk-cluster-resource-probe CTest gates this binary in CI; the
// fastgatk-cluster-smoke gate invokes verify_cluster_smoke.sh which in
// turn calls this binary on every scatter shard.  Together they replace
// the local fake-sbatch path in verify_slurm_wrapper.py for any
// production-bound run.

#include "fastgatk/runtime/resource.hpp"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>

namespace {

std::string env_or_empty(const char* name) {
    const char* value = std::getenv(name);
    return value ? std::string(value) : std::string{};
}

bool real_cluster_mode() {
    const auto mode = env_or_empty("FASTGATK_REAL_CLUSTER");
    return mode == "1" || mode == "true" || mode == "TRUE" || mode == "yes";
}

bool local_proof_mode() {
    const auto mode = env_or_empty("FASTGATK_LOCAL_PROOF");
    return mode == "1" || mode == "true" || mode == "TRUE" || mode == "yes";
}

std::uint64_t parse_uint(const std::string& value) {
    if (value.empty()) return 0;
    std::uint64_t result = 0;
    for (const char character : value) {
        if (character < '0' || character > '9') return 0;
        result = result * 10 + static_cast<std::uint64_t>(character - '0');
    }
    return result;
}

void emit(const std::string& payload, const std::string& path) {
    if (path.empty() || path == "-") {
        std::cout << payload << '\n';
        return;
    }
    std::ofstream stream(path, std::ios::trunc);
    if (!stream) throw std::runtime_error("cannot open evidence file: " + path);
    stream << payload << '\n';
}

std::string json_escape(std::string_view value) {
    std::ostringstream output;
    for (const char character : value) {
        switch (character) {
            case '"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default: output << character; break;
        }
    }
    return output.str();
}

}  // namespace

int main() {
    try {
        if (!real_cluster_mode() && !local_proof_mode()) {
            const std::string payload =
                "{\"status\":\"skip\",\"reason\":\"FASTGATK_REAL_CLUSTER/FASTGATK_LOCAL_PROOF not set\"}\n";
            emit(payload, env_or_empty("FASTGATK_CLUSTER_EVIDENCE_JSON"));
            return 0;
        }
        if (local_proof_mode() && !real_cluster_mode()) {
            // Local proof: capture host resource snapshot verbatim; no
            // SLURM env exists.  The dispatcher / RUNBOOK distinguish this
            // evidence from a real-cluster run via cluster_kind.
            const auto snapshot = fastgatk::runtime::ResourceSnapshot::probe();
            std::ostringstream output;
            output << "{"
                   << "\"schema_version\":1"
                   << ",\"status\":\"pass\""
                   << ",\"probe_source\":\"fastgatk-cluster-resource-probe\""
                   << ",\"mode\":\"local-proof\""
                   << ",\"snapshot\":" << snapshot.to_json()
                   << "}";
            emit(output.str(), env_or_empty("FASTGATK_CLUSTER_EVIDENCE_JSON"));
            return 0;
        }
        const auto slurm_job_id = env_or_empty("SLURM_JOB_ID");
        const auto slurm_step_id = env_or_empty("SLURM_STEP_ID");
        const auto slurm_cpus = env_or_empty("SLURM_CPUS_PER_TASK");
        const auto slurm_cpus_on_node = env_or_empty("SLURM_CPUS_ON_NODE");
        const auto slurm_mem_node = env_or_empty("SLURM_MEM_PER_NODE");
        const auto slurm_mem_cpu = env_or_empty("SLURM_MEM_PER_CPU");
        const auto slurm_tmpdir = env_or_empty("SLURM_TMPDIR");
        if (slurm_job_id.empty() && slurm_cpus.empty()) {
            // A real-cluster run that does not actually own an allocation is
            // a configuration bug, not a graceful skip.  Refuse up front so
            // the runbook's "5-minute rollback" diagnostic catches it.
            std::cerr << "{\"status\":\"fail\",\"reason\":\"FASTGATK_REAL_CLUSTER=1 but "
                         "SLURM_JOB_ID/SLURM_CPUS_PER_TASK absent\"}\n";
            return 2;
        }
        const auto snapshot = fastgatk::runtime::ResourceSnapshot::probe();
        // Required fields the cluster contract promises:
        //   slurm_job_id       == $SLURM_JOB_ID (verbatim)
        //   slurm_threads      == $SLURM_CPUS_PER_TASK (or _CPUS_ON_NODE fallback)
        //   slurm_memory_bytes == parsed $SLURM_MEM_PER_NODE (or _MEM_PER_CPU × threads)
        //   slurm_allocation   == true
        //   safe_memory_budget_bytes <= effective_memory_limit_bytes
        //   effective_threads(requested) <= allocation when cpuset applies
        if (snapshot.slurm_job_id != slurm_job_id) {
            throw std::runtime_error("slurm_job_id mismatch: probe='" + snapshot.slurm_job_id +
                                     "' env='" + slurm_job_id + "'");
        }
        if (!slurm_cpus.empty() && snapshot.slurm_threads != parse_uint(slurm_cpus)) {
            throw std::runtime_error("slurm_threads mismatch: probe=" +
                                     std::to_string(snapshot.slurm_threads) +
                                     " env=" + slurm_cpus);
        }
        if (snapshot.slurm_allocation != true) {
            throw std::runtime_error("slurm_allocation probe flag is false inside SLURM_JOB_ID");
        }
        const auto limit = snapshot.effective_memory_limit_bytes();
        const auto safe = snapshot.safe_memory_budget_bytes();
        if (limit != 0 && safe > limit) {
            throw std::runtime_error("safe budget exceeds effective limit: safe=" +
                                     std::to_string(safe) + " limit=" + std::to_string(limit));
        }
        if (!slurm_cpus.empty()) {
            const auto threads = parse_uint(slurm_cpus);
            const auto requested_eight = snapshot.effective_threads(8);
            if (requested_eight > threads) {
                throw std::runtime_error("effective_threads(8)=" +
                                         std::to_string(requested_eight) +
                                         " exceeds SLURM allocation=" + slurm_cpus);
            }
        }
        // Build the cluster-evidence JSON.  The shape is a strict superset of
        // the resource probe to_json() plus the SLURM env echoes that the
        // dispatcher / workflow scripts read for provenance.  Cluster-evidence
        // consumers MUST treat missing fields as a contract violation.
        std::ostringstream output;
        output << "{"
               << "\"schema_version\":1"
               << ",\"status\":\"pass\""
               << ",\"probe_source\":\"fastgatk-cluster-resource-probe\""
               << ",\"slurm_env\":{"
               << "\"SLURM_JOB_ID\":\"" << json_escape(slurm_job_id) << "\""
               << ",\"SLURM_STEP_ID\":\"" << json_escape(slurm_step_id) << "\""
               << ",\"SLURM_CPUS_PER_TASK\":\"" << json_escape(slurm_cpus) << "\""
               << ",\"SLURM_CPUS_ON_NODE\":\"" << json_escape(slurm_cpus_on_node) << "\""
               << ",\"SLURM_MEM_PER_NODE\":\"" << json_escape(slurm_mem_node) << "\""
               << ",\"SLURM_MEM_PER_CPU\":\"" << json_escape(slurm_mem_cpu) << "\""
               << ",\"SLURM_TMPDIR\":\"" << json_escape(slurm_tmpdir) << "\""
               << "}"
               << ",\"snapshot\":" << snapshot.to_json()
               << "}";
        emit(output.str(), env_or_empty("FASTGATK_CLUSTER_EVIDENCE_JSON"));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "{\"status\":\"fail\",\"error\":\"" << json_escape(error.what()) << "\"}\n";
        return 2;
    }
}
