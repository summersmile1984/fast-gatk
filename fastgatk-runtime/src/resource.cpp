#include "fastgatk/runtime/resource.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <thread>

#if defined(__linux__)
#include <sched.h>
#include <sys/statvfs.h>
#endif

namespace fastgatk::runtime {
namespace {

std::string env(const char* name) {
    const char* value = std::getenv(name);
    return value ? std::string(value) : std::string{};
}

std::uint64_t parse_uint(const std::string& value) {
    if (value.empty() || value == "max") return 0;
    std::uint64_t result = 0;
    const auto* first = value.data();
    const auto* last = first + value.size();
    const auto parsed = std::from_chars(first, last, result);
    return parsed.ec == std::errc{} ? result : 0;
}

bool parse_bool(const std::string& value) {
    return value == "1" || value == "true" || value == "TRUE" ||
           value == "yes" || value == "YES";
}

std::size_t visible_device_count(const std::string& value) {
    if (value.empty() || value == "NoDevFiles" || value == "void" ||
        value == "none" || value == "None")
        return 0;
    std::size_t count = 1;
    for (const char character : value) {
        if (character == ',') ++count;
    }
    return count;
}

std::pair<std::size_t, std::string> launcher_visible_devices() {
    // Prefer backend-specific selectors. SLURM_JOB_GPUS is intentionally last
    // because a launcher may also provide a narrower runtime selector. These
    // fields establish allocation provenance only; they do not create a
    // device context.
    constexpr std::array<std::pair<const char*, const char*>, 5> selectors{{
        {"CUDA_VISIBLE_DEVICES", "CUDA_VISIBLE_DEVICES"},
        {"HIP_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES"},
        {"ROCR_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES"},
        {"ZE_AFFINITY_MASK", "ZE_AFFINITY_MASK"},
        {"SLURM_JOB_GPUS", "SLURM_JOB_GPUS"},
    }};
    for (const auto& [name, label] : selectors) {
        const auto value = env(name);
        if (!value.empty()) return {visible_device_count(value), label};
    }
    return {0, {}};
}

std::uint64_t parse_slurm_memory(const std::string& value) {
    if (value.empty()) return 0;
    std::string number = value;
    std::uint64_t multiplier = 1024ULL * 1024ULL;  // SLURM memory values default to MB.
    const char suffix = number.back();
    if (suffix == 'K' || suffix == 'k') { multiplier = 1024ULL; number.pop_back(); }
    else if (suffix == 'G' || suffix == 'g') { multiplier = 1024ULL * 1024ULL * 1024ULL; number.pop_back(); }
    else if (suffix == 'T' || suffix == 't') { multiplier = 1024ULL * 1024ULL * 1024ULL * 1024ULL; number.pop_back(); }
    return parse_uint(number) * multiplier;
}

std::uint64_t read_limit_file(const std::filesystem::path& path) {
    std::ifstream input(path);
    std::string value;
    if (!input || !(input >> value)) return 0;
    return parse_uint(value);
}

std::filesystem::path cgroup_v2_path() {
    std::ifstream input("/proc/self/cgroup");
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("0::", 0) != 0) continue;
        std::filesystem::path relative = line.substr(3);
        // /proc/self/cgroup paths are absolute within the hierarchy; remove
        // the leading slash before appending to the mounted cgroup root.
        if (relative.is_absolute()) relative = relative.relative_path();
        return std::filesystem::path("/sys/fs/cgroup") / relative;
    }
    return std::filesystem::path("/sys/fs/cgroup");
}

std::uint64_t read_cgroup_value(const std::filesystem::path& root,
                                const char* filename) {
    const auto nested = read_limit_file(root / filename);
    if (nested != 0) return nested;
    if (root != std::filesystem::path("/sys/fs/cgroup"))
        return read_limit_file(std::filesystem::path("/sys/fs/cgroup") / filename);
    return 0;
}

std::size_t read_file_descriptor_limit() {
#if defined(__linux__)
    std::ifstream input("/proc/self/limits");
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("Max open files", 0) != 0) continue;
        std::istringstream fields(line);
        std::string word;
        while (fields >> word) {
            const auto value = parse_uint(word);
            if (value != 0) return static_cast<std::size_t>(value);
        }
        break;
    }
#endif
    return 0;
}

std::string json_escape(const std::string& value) {
    std::ostringstream output;
    for (const char character : value) {
        if (character == '\\' || character == '"') output << '\\';
        if (character == '\n') output << "\\n";
        else if (character == '\r') output << "\\r";
        else if (character == '\t') output << "\\t";
        else output << character;
    }
    return output.str();
}

std::uint64_t fraction(std::uint64_t value, std::uint64_t numerator,
                       std::uint64_t denominator) {
    if (value == 0) return 0;
    if (value > std::numeric_limits<std::uint64_t>::max() / numerator)
        return std::numeric_limits<std::uint64_t>::max() / denominator * numerator;
    return (value * numerator) / denominator;
}

std::uint32_t scaled_reads(std::uint64_t target, std::uint64_t bytes_per_read,
                           std::uint32_t requested) {
    if (target == 0 || bytes_per_read == 0)
        return requested == 0 ? std::numeric_limits<std::uint32_t>::max() : requested;
    const auto count = target / bytes_per_read;
    if (count == 0) return 1;
    const auto bounded = requested == 0 ? count : std::min<std::uint64_t>(count, requested);
    return static_cast<std::uint32_t>(std::min<std::uint64_t>(
        bounded, std::numeric_limits<std::uint32_t>::max()));
}

bool over_eighty(std::uint64_t used, std::uint64_t limit) {
    return limit != 0 && used > fraction(limit, 4, 5);
}

}  // namespace

BatchLimits AdaptiveController::initial(const ResourceBudget& budget,
                                        const WorkEstimate& work) const {
    const auto host_target = fraction(budget.host_target_bytes, 3, 5);
    const auto device_target = fraction(budget.device_free_bytes, 3, 5);
    auto inflight_target = host_target != 0 ? host_target : device_target;
    const auto scratch_target = fraction(budget.scratch_free_bytes, 4, 5);
    if (scratch_target != 0)
        inflight_target = inflight_target == 0 ? scratch_target :
            std::min(inflight_target, scratch_target);
    BatchLimits limits;
    limits.max_reads = std::max<std::uint32_t>(1, scaled_reads(
        host_target, work.host_bytes_per_read, work.reads));
    if (device_target != 0 && work.device_bytes_per_read != 0)
        limits.max_reads = std::min(limits.max_reads, scaled_reads(
            device_target, work.device_bytes_per_read, work.reads));
    if (inflight_target != 0 && work.inflight_bytes_per_read != 0)
        limits.max_reads = std::min(limits.max_reads, scaled_reads(
            inflight_target, work.inflight_bytes_per_read, work.reads));
    limits.max_haplotypes = std::max<std::uint32_t>(1, work.haplotypes);
    limits.max_host_bytes = host_target;
    limits.max_device_bytes = device_target;
    limits.max_inflight_bytes = inflight_target;
    return limits;
}

BatchLimits AdaptiveController::next(const RuntimeTelemetry& telemetry,
                                     const BatchLimits& current,
                                     Pressure pressure) const {
    BatchLimits next = current;
    const bool host_pressure = pressure == Pressure::HostMemory ||
        over_eighty(telemetry.host_bytes, current.max_host_bytes);
    const bool device_pressure = pressure == Pressure::DeviceMemory ||
        over_eighty(telemetry.device_bytes, current.max_device_bytes);
    const bool inflight_pressure = pressure == Pressure::IoBound ||
        over_eighty(telemetry.inflight_bytes, current.max_inflight_bytes);
    const bool scratch_pressure = pressure == Pressure::Scratch;
    const bool fd_pressure = pressure == Pressure::FileDescriptors;
    if (host_pressure || device_pressure || inflight_pressure || scratch_pressure || fd_pressure) {
        next.max_reads = std::max<std::uint32_t>(1, current.max_reads / 2);
        next.max_haplotypes = std::max<std::uint32_t>(1, current.max_haplotypes / 2);
        if (current.max_host_bytes != 0) next.max_host_bytes = std::max<std::uint64_t>(1, current.max_host_bytes / 2);
        if (current.max_device_bytes != 0) next.max_device_bytes = std::max<std::uint64_t>(1, current.max_device_bytes / 2);
        if (current.max_inflight_bytes != 0) next.max_inflight_bytes = std::max<std::uint64_t>(1, current.max_inflight_bytes / 2);
        return next;
    }
    if (pressure == Pressure::Normal && telemetry.compute_queue_empty &&
        telemetry.host_bytes < fraction(current.max_host_bytes, 1, 2) &&
        telemetry.device_bytes < fraction(current.max_device_bytes, 1, 2) &&
        telemetry.inflight_bytes < fraction(current.max_inflight_bytes, 1, 2)) {
        const auto grow = [](std::uint32_t value) {
            const auto delta = std::max<std::uint32_t>(1, value / 5);
            return value > std::numeric_limits<std::uint32_t>::max() - delta
                ? std::numeric_limits<std::uint32_t>::max() : value + delta;
        };
        next.max_reads = grow(current.max_reads);
        next.max_haplotypes = grow(current.max_haplotypes);
        const auto grow_bytes = [](std::uint64_t value) {
            const auto delta = std::max<std::uint64_t>(1, value / 5);
            return value > std::numeric_limits<std::uint64_t>::max() - delta
                ? std::numeric_limits<std::uint64_t>::max() : value + delta;
        };
        if (current.max_host_bytes != 0) next.max_host_bytes = grow_bytes(current.max_host_bytes);
        if (current.max_device_bytes != 0) next.max_device_bytes = grow_bytes(current.max_device_bytes);
        if (current.max_inflight_bytes != 0) next.max_inflight_bytes = grow_bytes(current.max_inflight_bytes);
    }
    return next;
}

bool AdaptiveController::should_spill(const RuntimeTelemetry& telemetry,
                                      const ResourceBudget& budget) const {
    return over_eighty(telemetry.scratch_bytes, budget.scratch_free_bytes) ||
           over_eighty(telemetry.scratch_inodes, budget.scratch_free_inodes);
}

bool AdaptiveController::should_pause_reader(const RuntimeTelemetry& telemetry,
                                             const ResourceBudget& budget) const {
    return over_eighty(telemetry.host_bytes, budget.host_target_bytes) ||
           over_eighty(telemetry.device_bytes, budget.device_free_bytes) ||
           over_eighty(telemetry.inflight_bytes, budget.host_target_bytes) ||
           over_eighty(telemetry.scratch_bytes, budget.scratch_free_bytes) ||
           over_eighty(telemetry.scratch_inodes, budget.scratch_free_inodes) ||
           (budget.file_descriptor_limit != 0 &&
            telemetry.open_file_descriptors > budget.file_descriptor_limit * 4 / 5);
}

ResourceSnapshot ResourceSnapshot::probe() {
    ResourceSnapshot snapshot;
#if defined(FASTGATK_RUNTIME_EXECUTION_SPACE)
    snapshot.backend = FASTGATK_RUNTIME_EXECUTION_SPACE;
#endif
#if defined(FASTGATK_RUNTIME_DEVICE_BACKEND) && FASTGATK_RUNTIME_DEVICE_BACKEND
    snapshot.device_backend_compiled = true;
#endif
    snapshot.hardware_threads = std::max<std::size_t>(1, std::thread::hardware_concurrency());
#if defined(__linux__)
    // libgomp may bind the initial thread to its first OpenMP place before
    // Kokkos is initialized.  That is a placement policy, not a scheduler
    // allocation: treating the one-place affinity as the machine allocation
    // silently turns GATK's default native PairHMM width (four) into one.
    // A SLURM allocation still supplies the authoritative cap below.  Only
    // ignore the transient OpenMP binding; ordinary cpuset affinity remains
    // a hard local limit.
    const auto omp_proc_bind = env("OMP_PROC_BIND");
    const bool openmp_binds_initial_thread = !omp_proc_bind.empty() &&
        omp_proc_bind != "0" && omp_proc_bind != "false" && omp_proc_bind != "FALSE";
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    if (!openmp_binds_initial_thread &&
        sched_getaffinity(0, sizeof(affinity), &affinity) == 0) {
        const auto visible = static_cast<std::size_t>(CPU_COUNT(&affinity));
        if (visible != 0) snapshot.hardware_threads = std::min(snapshot.hardware_threads, visible);
    }
#endif
    snapshot.slurm_job_id = env("SLURM_JOB_ID");
    snapshot.slurm_allocation = !snapshot.slurm_job_id.empty() || !env("SLURM_CPUS_PER_TASK").empty();
    snapshot.slurm_threads = static_cast<std::size_t>(parse_uint(env("SLURM_CPUS_PER_TASK")));
    if (snapshot.slurm_threads == 0) snapshot.slurm_threads = static_cast<std::size_t>(parse_uint(env("SLURM_CPUS_ON_NODE")));
    snapshot.slurm_memory_bytes = parse_slurm_memory(env("SLURM_MEM_PER_NODE"));
    if (snapshot.slurm_memory_bytes == 0) {
        const auto per_cpu = parse_slurm_memory(env("SLURM_MEM_PER_CPU"));
        if (per_cpu != 0 && snapshot.slurm_threads != 0) snapshot.slurm_memory_bytes = per_cpu * snapshot.slurm_threads;
    }
    snapshot.scratch_directory = env("SLURM_TMPDIR");
    if (snapshot.scratch_directory.empty()) snapshot.scratch_directory = env("TMPDIR");
    // Site-specific launchers may identify local NVMe and object-store inputs
    // explicitly.  Do not infer either flag from path spelling.
    snapshot.local_ssd = parse_bool(env("FASTGATK_LOCAL_SSD"));
    snapshot.remote_input = parse_bool(env("FASTGATK_REMOTE_INPUT"));
    snapshot.file_descriptor_limit = read_file_descriptor_limit();
#if defined(__linux__)
    if (!snapshot.scratch_directory.empty()) {
        struct statvfs filesystem_stats{};
        if (statvfs(snapshot.scratch_directory.c_str(), &filesystem_stats) == 0) {
            if (filesystem_stats.f_bavail != 0 && filesystem_stats.f_frsize <=
                std::numeric_limits<std::uint64_t>::max() / filesystem_stats.f_bavail)
                snapshot.scratch_free_bytes =
                    static_cast<std::uint64_t>(filesystem_stats.f_bavail) * filesystem_stats.f_frsize;
            snapshot.scratch_free_inodes = filesystem_stats.f_favail;
        }
    }
#endif

    const auto cgroup_root = cgroup_v2_path();
    snapshot.memory_limit_bytes = read_cgroup_value(cgroup_root, "memory.max");
    snapshot.memory_current_bytes = read_cgroup_value(cgroup_root, "memory.current");
    snapshot.cgroup_memory_limited = snapshot.memory_limit_bytes != 0;

    // Kokkos deliberately has no portable free-device-memory API. Preserve
    // launcher values as provenance, but admit them to the controller only
    // when this binary has a device execution space and an assignment/memory
    // declaration is visible. A host build with CUDA_VISIBLE_DEVICES is common
    // on shared nodes and must fail closed.
    snapshot.reported_device_free_bytes = parse_uint(env("FASTGATK_DEVICE_FREE_BYTES"));
    snapshot.reported_device_memory_bytes = parse_uint(env("FASTGATK_DEVICE_MEMORY_BYTES"));
    if (snapshot.reported_device_memory_bytes == 0)
        snapshot.reported_device_memory_bytes = parse_uint(env("FASTGATK_DEVICE_TOTAL_BYTES"));
    const auto [visible_devices, visibility_source] = launcher_visible_devices();
    snapshot.visible_gpus = visible_devices;
    snapshot.device_visibility_source = visibility_source;
    const bool device_reported = snapshot.visible_gpus != 0 ||
        snapshot.reported_device_memory_bytes != 0 || snapshot.reported_device_free_bytes != 0;
    snapshot.device_assignment_available = snapshot.device_backend_compiled && device_reported;
    if (snapshot.device_assignment_available) {
        snapshot.device_memory_bytes = snapshot.reported_device_memory_bytes;
        snapshot.device_free_bytes = snapshot.reported_device_free_bytes;
    } else if (device_reported) {
        snapshot.device_telemetry_rejected = true;
    }
    return snapshot;
}

std::size_t ResourceSnapshot::effective_threads(std::size_t requested) const {
    const std::size_t allocation = slurm_threads == 0 ? hardware_threads : slurm_threads;
    if (requested == 0) return std::max<std::size_t>(1, allocation);
    return std::max<std::size_t>(1, std::min(requested, allocation));
}

std::uint64_t ResourceSnapshot::effective_memory_limit_bytes() const {
    if (memory_limit_bytes == 0) return slurm_memory_bytes;
    if (slurm_memory_bytes == 0) return memory_limit_bytes;
    return std::min(memory_limit_bytes, slurm_memory_bytes);
}

std::uint64_t ResourceSnapshot::safe_memory_budget_bytes() const {
    const auto limit = effective_memory_limit_bytes();
    if (limit == 0) return 0;
    // Keep arithmetic overflow-safe for synthetic scheduler limits.
    return limit > std::numeric_limits<std::uint64_t>::max() / 8
        ? std::numeric_limits<std::uint64_t>::max() / 10 * 8
        : (limit * 8) / 10;
}

bool ResourceSnapshot::memory_budget_allows(std::uint64_t owned_bytes) const {
    const auto budget = safe_memory_budget_bytes();
    return budget == 0 || owned_bytes <= budget;
}

ResourceBudget ResourceSnapshot::budget() const {
    ResourceBudget result;
    result.host_target_bytes = safe_memory_budget_bytes();
    result.host_hard_bytes = effective_memory_limit_bytes();
    result.device_free_bytes = device_free_bytes;
    result.device_hard_bytes = device_memory_bytes;
    const auto directory = scratch_directory.empty() ? std::filesystem::path{} :
                           std::filesystem::path(scratch_directory);
    if (!directory.empty()) {
        std::error_code error;
        const auto space = std::filesystem::space(directory, error);
        if (!error) result.scratch_free_bytes = space.available;
    }
    if (result.scratch_free_bytes == 0) result.scratch_free_bytes = scratch_free_bytes;
    result.scratch_hard_bytes = result.scratch_free_bytes;
    result.scratch_free_inodes = scratch_free_inodes;
    result.file_descriptor_limit = file_descriptor_limit;
    result.cpu_threads = effective_threads();
    result.local_ssd = local_ssd;
    result.remote_input = remote_input;
    return result;
}

std::string ResourceSnapshot::to_json() const {
    std::ostringstream output;
    output << "{\"backend\":\"" << json_escape(backend)
           << "\",\"execution_space\":\"" << json_escape(backend)
           << "\",\"device_backend_compiled\":" << (device_backend_compiled ? "true" : "false")
           << ",\"device_assignment_available\":" << (device_assignment_available ? "true" : "false")
           << ",\"device_telemetry_rejected\":" << (device_telemetry_rejected ? "true" : "false")
           << ",\"device_visibility_source\":\"" << json_escape(device_visibility_source)
           << "\",\"device_runtime_probe\":\"not-run (resource probe never initializes a device context)\""
           << ",\"memory_limit_bytes\":" << effective_memory_limit_bytes()
           << ",\"safe_memory_budget_bytes\":" << safe_memory_budget_bytes()
           << ",\"scratch_free_bytes\":" << budget().scratch_free_bytes
           << ",\"host_hard_bytes\":" << budget().host_hard_bytes
           << ",\"device_memory_bytes\":" << device_memory_bytes
           << ",\"device_free_bytes\":" << device_free_bytes
           << ",\"reported_device_memory_bytes\":" << reported_device_memory_bytes
           << ",\"reported_device_free_bytes\":" << reported_device_free_bytes
           << ",\"scratch_hard_bytes\":" << budget().scratch_hard_bytes
           << ",\"scratch_free_inodes\":" << budget().scratch_free_inodes
           << ",\"file_descriptor_limit\":" << budget().file_descriptor_limit
           << ",\"memory_current_bytes\":" << memory_current_bytes
           << ",\"hardware_threads\":" << hardware_threads
           << ",\"slurm_threads\":" << slurm_threads
           << ",\"visible_gpus\":" << visible_gpus
           << ",\"cpu_threads\":" << budget().cpu_threads
           << ",\"local_ssd\":" << (local_ssd ? "true" : "false")
           << ",\"remote_input\":" << (remote_input ? "true" : "false")
           << ",\"cgroup_memory_limited\":" << (cgroup_memory_limited ? "true" : "false")
           << ",\"slurm_allocation\":" << (slurm_allocation ? "true" : "false")
           << ",\"slurm_job_id\":\"" << json_escape(slurm_job_id)
           << "\",\"scratch_directory\":\"" << json_escape(scratch_directory) << "\"}";
    return output.str();
}

}  // namespace fastgatk::runtime
