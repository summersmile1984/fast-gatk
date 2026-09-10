#pragma once

#include <cstddef>
#include <cstdint>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <optional>
#include <string>
#include <utility>

namespace fastgatk::runtime {

// The runtime-facing budget is intentionally independent of Kokkos and
// HTSlib. A zero field means that the corresponding limit is unknown and the
// controller must not manufacture a host-wide limit from it.
struct ResourceBudget {
    std::uint64_t host_target_bytes = 0;
    std::uint64_t device_free_bytes = 0;
    std::uint64_t scratch_free_bytes = 0;
    std::uint64_t scratch_free_inodes = 0;
    std::size_t file_descriptor_limit = 0;
    // Hard limits and scheduler/provenance flags are additive to the original
    // five-field aggregate so existing callers remain source compatible.
    // A zero hard limit means that the boundary is unknown, never unlimited.
    std::uint64_t host_hard_bytes = 0;
    std::uint64_t device_hard_bytes = 0;
    std::uint64_t scratch_hard_bytes = 0;
    std::size_t cpu_threads = 0;
    bool local_ssd = false;
    bool remote_input = false;
};

struct WorkEstimate {
    std::uint32_t reads = 1;
    std::uint32_t haplotypes = 1;
    std::uint64_t host_bytes_per_read = 0;
    std::uint64_t device_bytes_per_read = 0;
    std::uint64_t inflight_bytes_per_read = 0;
};

struct BatchLimits {
    std::uint32_t max_reads = 1;
    std::uint32_t max_haplotypes = 1;
    std::uint64_t max_host_bytes = 0;
    std::uint64_t max_device_bytes = 0;
    std::uint64_t max_inflight_bytes = 0;
};

enum class Pressure {
    Normal,
    HostMemory,
    DeviceMemory,
    Scratch,
    FileDescriptors,
    IoBound,
};

struct RuntimeTelemetry {
    std::uint64_t host_bytes = 0;
    std::uint64_t device_bytes = 0;
    std::uint64_t inflight_bytes = 0;
    std::uint64_t scratch_bytes = 0;
    std::uint64_t scratch_inodes = 0;
    std::size_t open_file_descriptors = 0;
    bool compute_queue_empty = false;
    bool reader_paused = false;
};

// Conservative controller used at safe batch boundaries. It never increases
// a limit under pressure and never allows a zero-byte estimate to turn into an
// accidental divide-by-zero or unbounded allocation.
class AdaptiveController {
public:
    BatchLimits initial(const ResourceBudget& budget, const WorkEstimate& work) const;
    BatchLimits next(const RuntimeTelemetry& telemetry,
                     const BatchLimits& current,
                     Pressure pressure) const;
    bool should_spill(const RuntimeTelemetry& telemetry,
                      const ResourceBudget& budget) const;
    bool should_pause_reader(const RuntimeTelemetry& telemetry,
                             const ResourceBudget& budget) const;
};

// A byte-capacity queue used between decode, compute and encode stages. The
// queue accounts bytes rather than item count, so a single large read/chunk
// cannot bypass backpressure. `push` blocks until capacity is available and
// returns false after close; `pop` returns nullopt after close and drain.
template<class Record>
class BoundedByteQueue {
public:
    explicit BoundedByteQueue(std::uint64_t capacity_bytes)
        : capacity_bytes_(capacity_bytes) {}

    BoundedByteQueue(const BoundedByteQueue&) = delete;
    BoundedByteQueue& operator=(const BoundedByteQueue&) = delete;

    bool push(Record record, std::uint64_t bytes) {
        if (bytes > capacity_bytes_) return false;
        std::unique_lock lock(mutex_);
        not_full_.wait(lock, [&] {
            return closed_ || queued_bytes_ <= capacity_bytes_ - bytes;
        });
        if (closed_) return false;
        queued_bytes_ += bytes;
        if (queued_bytes_ > peak_queued_bytes_) peak_queued_bytes_ = queued_bytes_;
        queue_.emplace_back(std::move(record), bytes);
        not_empty_.notify_one();
        return true;
    }

    std::optional<Record> pop(std::uint64_t* bytes = nullptr) {
        std::unique_lock lock(mutex_);
        not_empty_.wait(lock, [&] { return closed_ || !queue_.empty(); });
        if (queue_.empty()) return std::nullopt;
        auto entry = std::move(queue_.front());
        queue_.pop_front();
        queued_bytes_ -= entry.second;
        if (bytes != nullptr) *bytes = entry.second;
        not_full_.notify_all();
        return std::move(entry.first);
    }

    void close() {
        std::lock_guard lock(mutex_);
        closed_ = true;
        not_full_.notify_all();
        not_empty_.notify_all();
    }

    bool closed() const {
        std::lock_guard lock(mutex_);
        return closed_;
    }

    std::uint64_t queued_bytes() const {
        std::lock_guard lock(mutex_);
        return queued_bytes_;
    }

    // Highest occupancy observed since construction.  This is maintained
    // while holding the queue mutex, so a fast consumer cannot erase the
    // evidence of a producer-side peak before telemetry samples it.
    std::uint64_t peak_queued_bytes() const {
        std::lock_guard lock(mutex_);
        return peak_queued_bytes_;
    }

    std::size_t size() const {
        std::lock_guard lock(mutex_);
        return queue_.size();
    }

private:
    const std::uint64_t capacity_bytes_;
    mutable std::mutex mutex_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    std::deque<std::pair<Record, std::uint64_t>> queue_;
    std::uint64_t queued_bytes_ = 0;
    std::uint64_t peak_queued_bytes_ = 0;
    bool closed_ = false;
};

// A read-only view of the allocation visible to the native process.  It is
// deliberately independent of Kokkos so the dispatcher and fallback path can
// probe resources before initializing an execution space.
struct ResourceSnapshot {
    std::uint64_t memory_limit_bytes = 0;
    std::uint64_t memory_current_bytes = 0;
    std::uint64_t slurm_memory_bytes = 0;
    // Optional accelerator memory that is usable by *this binary*.  Kokkos
    // itself does not expose a portable free-memory query, so a missing value
    // remains zero and the controller must use its own bounded policy.  A
    // launcher-reported device budget is deliberately not admitted here for a
    // host-only binary; see reported_device_* below for provenance.
    std::uint64_t device_memory_bytes = 0;
    std::uint64_t device_free_bytes = 0;
    std::uint64_t reported_device_memory_bytes = 0;
    std::uint64_t reported_device_free_bytes = 0;
    std::uint64_t scratch_free_bytes = 0;
    std::uint64_t scratch_free_inodes = 0;
    std::size_t file_descriptor_limit = 0;
    std::size_t hardware_threads = 1;
    std::size_t slurm_threads = 0;
    std::size_t visible_gpus = 0;
    // Compile provenance and scheduler assignment are separate facts.  In
    // particular, CUDA_VISIBLE_DEVICES must not make an OpenMP/Serial binary
    // claim that it can allocate device Views.
    bool device_backend_compiled = false;
    bool device_assignment_available = false;
    bool device_telemetry_rejected = false;
    bool cgroup_memory_limited = false;
    bool slurm_allocation = false;
    bool local_ssd = false;
    bool remote_input = false;
    std::string slurm_job_id;
    std::string scratch_directory;
    std::string backend = "Host";
    std::string device_visibility_source;

    static ResourceSnapshot probe();

    // requested==0 means use the allocation.  The result never exceeds the
    // visible SLURM CPU allocation and is at least one.
    std::size_t effective_threads(std::size_t requested = 0) const;

    // Conservative hard limit used by batch controllers.  Zero means unknown.
    std::uint64_t effective_memory_limit_bytes() const;

    // Conservative 80% budget for owned Host/device staging. Zero means the
    // allocation limit is unknown and callers should use their own policy.
    std::uint64_t safe_memory_budget_bytes() const;

    bool memory_budget_allows(std::uint64_t owned_bytes) const;

    ResourceBudget budget() const;

    std::string to_json() const;
};

}  // namespace fastgatk::runtime
