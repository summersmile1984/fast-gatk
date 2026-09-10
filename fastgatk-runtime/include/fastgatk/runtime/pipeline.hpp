#pragma once

#include "fastgatk/runtime/resource.hpp"

#include <atomic>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace fastgatk::runtime {

// A bounded decode -> compute -> encode pipeline shared by native tools.
//
// The callbacks deliberately operate on value types owned by the queues.  A
// The Decode callback returns nullopt at end-of-input; Compute and Encode may
// return nullopt to drop an item (for example a read rejected by a Host filter)
// without changing ordering of items that reach the sink.  Byte functions are
// evaluated before an item enters the next queue, so every
// stage is subject to the same byte-capacity backpressure contract.  The
// implementation is backend agnostic: the compute callback may invoke a
// Kokkos execution space, while the pipeline itself only owns Host threads
// and synchronization.
template<class Decoded, class Computed, class Encoded>
class ThreeStagePipeline {
public:
    struct Limits {
        std::uint64_t decoded_bytes = 0;
        std::uint64_t computed_bytes = 0;
        std::uint64_t encoded_bytes = 0;
    };

    struct Metrics {
        std::uint64_t decoded_items = 0;
        std::uint64_t computed_items = 0;
        std::uint64_t encoded_items = 0;
        std::uint64_t decoded_bytes = 0;
        std::uint64_t computed_bytes = 0;
        std::uint64_t encoded_bytes = 0;
        std::uint64_t peak_decoded_bytes = 0;
        std::uint64_t peak_computed_bytes = 0;
        std::uint64_t peak_encoded_bytes = 0;
    };

    using Decode = std::function<std::optional<Decoded>()>;
    using Compute = std::function<std::optional<Computed>(Decoded)>;
    using Encode = std::function<std::optional<Encoded>(Computed)>;
    using Sink = std::function<void(Encoded)>;
    using ByteSize = std::function<std::uint64_t(const Decoded&)>;
    using ComputedByteSize = std::function<std::uint64_t(const Computed&)>;
    using EncodedByteSize = std::function<std::uint64_t(const Encoded&)>;

    ThreeStagePipeline(Limits limits,
                       Decode decode,
                       Compute compute,
                       Encode encode,
                       Sink sink,
                       ByteSize decoded_size,
                       ComputedByteSize computed_size,
                       EncodedByteSize encoded_size)
        : decoded_(require_capacity(limits.decoded_bytes, "decoded")),
          computed_(require_capacity(limits.computed_bytes, "computed")),
          encoded_(require_capacity(limits.encoded_bytes, "encoded")),
          decode_(std::move(decode)),
          compute_(std::move(compute)),
          encode_(std::move(encode)),
          sink_(std::move(sink)),
          decoded_size_(std::move(decoded_size)),
          computed_size_(std::move(computed_size)),
          encoded_size_(std::move(encoded_size)) {
        if (!decode_ || !compute_ || !encode_ || !sink_ ||
            !decoded_size_ || !computed_size_ || !encoded_size_)
            throw std::invalid_argument("ThreeStagePipeline callbacks must be non-empty");
    }

    ThreeStagePipeline(const ThreeStagePipeline&) = delete;
    ThreeStagePipeline& operator=(const ThreeStagePipeline&) = delete;

    Metrics run() {
        Metrics metrics;
        Metrics decoder_metrics;
        Metrics worker_metrics;
        Metrics encoder_metrics;
        std::exception_ptr failure;
        std::mutex failure_mutex;
        const auto fail = [&](std::exception_ptr error) {
            std::lock_guard lock(failure_mutex);
            if (!failure) failure = error;
            decoded_.close();
            computed_.close();
            encoded_.close();
        };
        std::thread decoder([&] {
            try {
                while (true) {
                    auto item = decode_();
                    if (!item.has_value()) break;
                    const auto bytes = decoded_size_(*item);
                    if (!decoded_.push(std::move(*item), bytes)) {
                        throw std::runtime_error(
                            "RESOURCE_EXHAUSTED: decoded pipeline queue is closed or item exceeds byte capacity");
                    }
                    ++decoder_metrics.decoded_items;
                    decoder_metrics.decoded_bytes += bytes;
                    decoder_metrics.peak_decoded_bytes = decoded_.peak_queued_bytes();
                }
                decoded_.close();
            } catch (...) {
                fail(std::current_exception());
            }
        });

        std::thread worker([&] {
            try {
                std::uint64_t bytes = 0;
                while (auto item = decoded_.pop(&bytes)) {
                    auto result = compute_(std::move(*item));
                    if (!result.has_value()) continue;
                    const auto result_bytes = computed_size_(*result);
                    if (!computed_.push(std::move(*result), result_bytes)) {
                        throw std::runtime_error(
                            "RESOURCE_EXHAUSTED: computed pipeline queue is closed or item exceeds byte capacity");
                    }
                    ++worker_metrics.computed_items;
                    worker_metrics.computed_bytes += result_bytes;
                    worker_metrics.peak_computed_bytes = computed_.peak_queued_bytes();
                }
                computed_.close();
            } catch (...) {
                fail(std::current_exception());
            }
        });

        std::thread encoder([&] {
            try {
                std::uint64_t bytes = 0;
                while (auto item = computed_.pop(&bytes)) {
                    auto result = encode_(std::move(*item));
                    if (!result.has_value()) continue;
                    const auto result_bytes = encoded_size_(*result);
                    if (!encoded_.push(std::move(*result), result_bytes)) {
                        throw std::runtime_error(
                            "RESOURCE_EXHAUSTED: encoded pipeline queue is closed or item exceeds byte capacity");
                    }
                    ++encoder_metrics.encoded_items;
                    encoder_metrics.encoded_bytes += result_bytes;
                    encoder_metrics.peak_encoded_bytes = encoded_.peak_queued_bytes();
                }
                encoded_.close();
            } catch (...) {
                fail(std::current_exception());
            }
        });

        // Drain the final queue on the caller thread.  This is intentionally
        // the only place the output sink runs, giving all tools a deterministic
        // single-writer boundary even when compute callbacks use Kokkos.
        try {
            std::uint64_t bytes = 0;
            while (auto item = encoded_.pop(&bytes)) sink_(std::move(*item));
        } catch (...) {
            fail(std::current_exception());
        }

        decoder.join();
        worker.join();
        encoder.join();
        if (failure) std::rethrow_exception(failure);
        metrics.decoded_items = decoder_metrics.decoded_items;
        metrics.computed_items = worker_metrics.computed_items;
        metrics.encoded_items = encoder_metrics.encoded_items;
        metrics.decoded_bytes = decoder_metrics.decoded_bytes;
        metrics.computed_bytes = worker_metrics.computed_bytes;
        metrics.encoded_bytes = encoder_metrics.encoded_bytes;
        metrics.peak_decoded_bytes = decoder_metrics.peak_decoded_bytes;
        metrics.peak_computed_bytes = worker_metrics.peak_computed_bytes;
        metrics.peak_encoded_bytes = encoder_metrics.peak_encoded_bytes;
        return metrics;
    }

private:
    static std::uint64_t require_capacity(std::uint64_t capacity, const char* stage) {
        if (capacity == 0)
            throw std::invalid_argument(std::string("ThreeStagePipeline ") + stage +
                                        " queue capacity must be positive");
        return capacity;
    }

    BoundedByteQueue<Decoded> decoded_;
    BoundedByteQueue<Computed> computed_;
    BoundedByteQueue<Encoded> encoded_;
    Decode decode_;
    Compute compute_;
    Encode encode_;
    Sink sink_;
    ByteSize decoded_size_;
    ComputedByteSize computed_size_;
    EncodedByteSize encoded_size_;
};

}  // namespace fastgatk::runtime
