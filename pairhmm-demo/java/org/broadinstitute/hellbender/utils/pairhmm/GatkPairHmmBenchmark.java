package org.broadinstitute.hellbender.utils.pairhmm;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.CyclicBarrier;

import org.broadinstitute.hellbender.utils.QualityUtils;
import com.intel.gkl.pairhmm.IntelPairHmm;
import com.intel.gkl.pairhmm.IntelPairHmmOMP;
import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;

/**
 * Calls the actual GATK 4.6.2.0 LoglessPairHMM class, not a reimplementation.
 * It deliberately uses the package-private testing entry point so the benchmark
 * measures the same one-pair kernel as the native demonstration.
 */
public final class GatkPairHmmBenchmark {
    private static final class SplitMix64 {
        long state;
        SplitMix64(long state) { this.state = state; }
        long next() {
            state += 0x9e3779b97f4a7c15L;
            long z = state;
            z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
            z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
            return z ^ (z >>> 31);
        }
        long mod(long n) { return Long.remainderUnsigned(next(), n); }
    }

    private static final class Pair {
        byte[] hap, read, readQual, ins, del, gcp;
    }

    private static int intArg(String[] args, String key, int fallback) {
        String prefix = key + "=";
        for (String arg : args) if (arg.startsWith(prefix)) return Integer.parseInt(arg.substring(prefix.length()));
        return fallback;
    }

    private static long longArg(String[] args, String key, long fallback) {
        String prefix = key + "=";
        for (String arg : args) if (arg.startsWith(prefix)) return Long.parseLong(arg.substring(prefix.length()));
        return fallback;
    }

    private static String stringArg(String[] args, String key, String fallback) {
        String prefix = key + "=";
        for (String arg : args) if (arg.startsWith(prefix)) return arg.substring(prefix.length());
        return fallback;
    }

    private static List<Pair> inputs(int count, int readLength, int hapLength, long seed) {
        final byte[] bases = {(byte)'A', (byte)'C', (byte)'G', (byte)'T'};
        final SplitMix64 rng = new SplitMix64(seed);
        final List<Pair> result = new ArrayList<>(count);
        for (int k = 0; k < count; ++k) {
            Pair p = new Pair();
            p.hap = new byte[hapLength]; p.read = new byte[readLength];
            p.readQual = new byte[readLength]; p.ins = new byte[readLength];
            p.del = new byte[readLength]; p.gcp = new byte[readLength];
            for (int i = 0; i < hapLength; ++i) p.hap[i] = bases[(int)rng.mod(4)];
            for (int i = 0; i < readLength; ++i) {
                p.read[i] = rng.mod(10) == 0 ? bases[(int)rng.mod(4)] : p.hap[i % hapLength];
                p.readQual[i] = (byte)(25 + rng.mod(16));
                p.ins[i] = (byte)(35 + rng.mod(8));
                p.del[i] = (byte)(35 + rng.mod(8));
                p.gcp[i] = (byte)(10 + rng.mod(10));
            }
            result.add(p);
        }
        return result;
    }

    public static void main(String[] args) throws Exception {
        final int count = intArg(args, "--pairs", 256);
        final int readLength = intArg(args, "--read-len", 150);
        final int hapLength = intArg(args, "--hap-len", 160);
        final int threads = intArg(args, "--threads", 1);
        final int iterations = intArg(args, "--iterations", 5);
        final long seed = longArg(args, "--seed", 42L);
        final String valuesOut = stringArg(args, "--values-out", "");
        final String sumsOut = stringArg(args, "--sums-out", "");
        final String tablesOut = stringArg(args, "--tables-out", "");
        final List<Pair> inputs = inputs(count, readLength, hapLength, seed);
        if (!tablesOut.isEmpty()) writeTables(tablesOut);
        final String mode = stringArg(args, "--mode", "java");
        final String workload = stringArg(args, "--workload", "independent");
        if (mode.equals("gatk-avx") || mode.equals("gatk-avx-omp")) {
            runNative(mode, inputs, threads, iterations, valuesOut);
            return;
        }
        if (workload.equals("matrix8")) {
            if (threads != 1) throw new IllegalArgumentException("matrix8 Java oracle currently requires --threads=1");
            runJavaMatrix8(inputs, readLength, hapLength, iterations, valuesOut, sumsOut);
            return;
        }
        if (!workload.equals("independent")) throw new IllegalArgumentException("unknown workload: " + workload);
        if (threads > 1) {
            runJavaParallel(inputs, readLength, hapLength, threads, iterations, valuesOut);
            return;
        }
        final LoglessPairHMM hmm = new LoglessPairHMM();
        hmm.initialize(readLength, hapLength);
        for (Pair p : inputs) hmm.computeReadLikelihoodGivenHaplotypeLog10(p.hap, p.read, p.readQual, p.ins, p.del, p.gcp, true, null);

        final long start = System.nanoTime();
        final double[] values = new double[count];
        for (int iteration = 0; iteration < iterations; ++iteration) {
            for (int i = 0; i < count; ++i) {
                Pair p = inputs.get(i);
                values[i] = hmm.computeReadLikelihoodGivenHaplotypeLog10(p.hap, p.read, p.readQual, p.ins, p.del, p.gcp, true, null);
            }
        }
        final double seconds = (System.nanoTime() - start) * 1.0e-9;
        double checksum = 0.0;
        for (double value : values) checksum += value;
        System.out.printf("{\"mode\":\"gatk-java-logless\",\"pairs\":%d,\"read_len\":%d,\"hap_len\":%d,\"threads\":%d,\"iterations\":%d,\"seconds\":%.12g,\"pairs_per_second\":%.12g,\"checksum\":%.12g}%n", count, readLength, hapLength, threads, iterations, seconds, count * (double)iterations / seconds, checksum);
        if (!valuesOut.isEmpty()) try (BufferedWriter writer = new BufferedWriter(new FileWriter(valuesOut))) {
            for (double value : values) writer.write(Long.toHexString(Double.doubleToRawLongBits(value)) + "\n");
        }
    }

    private static void runJavaMatrix8(List<Pair> inputs, int readLength, int hapLength,
                                       int iterations, String valuesOut, String sumsOut) throws IOException {
        final int block = Math.min(8, inputs.size());
        int pairCount = 0;
        for (int base = 0; base < inputs.size(); base += block) {
            final int n = Math.min(block, inputs.size() - base);
            pairCount += n * n;
        }
        final int[] readIds = new int[pairCount];
        final int[] hapIds = new int[pairCount];
        int cursor = 0;
        for (int base = 0; base < inputs.size(); base += block) {
            final int n = Math.min(block, inputs.size() - base);
            for (int read = 0; read < n; ++read) {
                for (int hap = 0; hap < n; ++hap) {
                    readIds[cursor] = base + read;
                    hapIds[cursor] = base + hap;
                    ++cursor;
                }
            }
        }
        final LoglessPairHMM hmm = new LoglessPairHMM();
        hmm.initialize(readLength, hapLength);
        final double[] values = new double[pairCount];
        final double[] scaledSums = new double[pairCount];
        for (int pair = 0; pair < pairCount; ++pair) {
            final Pair read = inputs.get(readIds[pair]);
            final Pair hap = inputs.get(hapIds[pair]);
            values[pair] = hmm.computeReadLikelihoodGivenHaplotypeLog10(
                    hap.hap, read.read, read.readQual, read.ins, read.del, read.gcp, true, null);
        }
        final long start = System.nanoTime();
        for (int iteration = 0; iteration < iterations; ++iteration) {
            for (int pair = 0; pair < pairCount; ++pair) {
                final Pair read = inputs.get(readIds[pair]);
                final Pair hap = inputs.get(hapIds[pair]);
                values[pair] = hmm.computeReadLikelihoodGivenHaplotypeLog10(
                        hap.hap, read.read, read.readQual, read.ins, read.del, read.gcp, true, null);
                double sum = 0.0;
                for (int j = 1; j <= hapLength; ++j)
                    sum += hmm.matchMatrix[readLength][j] + hmm.insertionMatrix[readLength][j];
                scaledSums[pair] = sum;
            }
        }
        final double seconds = (System.nanoTime() - start) * 1.0e-9;
        double checksum = 0.0;
        for (double value : values) checksum += value;
        System.out.printf("{\"mode\":\"gatk-java-logless\",\"workload\":\"matrix8\","
                        + "\"records\":%d,\"computed_pairs\":%d,\"threads\":1,\"iterations\":%d,"
                        + "\"seconds\":%.12g,\"pairs_per_second\":%.12g,\"checksum\":%.12g}%n",
                inputs.size(), pairCount, iterations, seconds,
                pairCount * (double)iterations / seconds, checksum);
        if (!valuesOut.isEmpty()) try (BufferedWriter writer = new BufferedWriter(new FileWriter(valuesOut))) {
            for (double value : values) writer.write(Long.toHexString(Double.doubleToRawLongBits(value)) + "\n");
        }
        if (!sumsOut.isEmpty()) try (BufferedWriter writer = new BufferedWriter(new FileWriter(sumsOut))) {
            for (double value : scaledSums) writer.write(Long.toHexString(Double.doubleToRawLongBits(value)) + "\n");
        }
    }

    private static void runNative(String mode, List<Pair> inputs, int threads, int iterations, String valuesOut) throws IOException {
        final PairHMMNativeArguments arguments = new PairHMMNativeArguments();
        arguments.maxNumberOfThreads = threads;
        arguments.useDoublePrecision = true;
        final PairHMMNativeBinding binding = mode.equals("gatk-avx-omp") ? new IntelPairHmmOMP() : new IntelPairHmm();
        if (!binding.load(null)) throw new IllegalStateException("GATK Intel GKL PairHMM could not load on this host");
        binding.initialize(arguments);
        final double[] values = new double[inputs.size()];
        final int block = Math.min(8, inputs.size());
        computeNativeBlocks(binding, inputs, values, block);
        final long start = System.nanoTime();
        for (int iteration = 0; iteration < iterations; ++iteration) computeNativeBlocks(binding, inputs, values, block);
        final double seconds = (System.nanoTime() - start) * 1.0e-9;
        double checksum = 0.0; for (double value : values) checksum += value;
        System.out.printf("{\"mode\":\"%s\",\"pairs\":%d,\"threads\":%d,\"iterations\":%d,\"seconds\":%.12g,\"pairs_per_second\":%.12g,\"checksum\":%.12g}%n", mode, inputs.size(), threads, iterations, seconds, inputs.size() * (double)iterations / seconds, checksum);
        if (!valuesOut.isEmpty()) try (BufferedWriter writer = new BufferedWriter(new FileWriter(valuesOut))) {
            for (double value : values) writer.write(Long.toHexString(Double.doubleToRawLongBits(value)) + "\n");
        }
        binding.done();
    }

    private static void runJavaParallel(List<Pair> inputs, int readLength, int hapLength, int threads, int iterations, String valuesOut) throws InterruptedException, IOException {
        final double[] values = new double[inputs.size()];
        final int actual = Math.min(threads, inputs.size());
        final AtomicInteger next = new AtomicInteger();
        final CyclicBarrier warmupBarrier = new CyclicBarrier(actual + 1);
        final CyclicBarrier startBarrier = new CyclicBarrier(actual + 1);
        final CyclicBarrier doneBarrier = new CyclicBarrier(actual + 1);
        final Thread[] workers = new Thread[actual];
        for (int t = 0; t < actual; ++t) {
            workers[t] = new Thread(() -> {
                final LoglessPairHMM hmm = new LoglessPairHMM();
                hmm.initialize(readLength, hapLength);
                try {
                    while (true) {
                        final int index = next.getAndIncrement();
                        if (index >= inputs.size()) break;
                        final Pair p = inputs.get(index);
                        values[index] = hmm.computeReadLikelihoodGivenHaplotypeLog10(
                                p.hap, p.read, p.readQual, p.ins, p.del, p.gcp, true, null);
                    }
                    warmupBarrier.await();
                    for (int iteration = 0; iteration < iterations; ++iteration) {
                        startBarrier.await();
                        while (true) {
                            final int index = next.getAndIncrement();
                            if (index >= inputs.size()) break;
                            final Pair p = inputs.get(index);
                            values[index] = hmm.computeReadLikelihoodGivenHaplotypeLog10(p.hap, p.read, p.readQual, p.ins, p.del, p.gcp, true, null);
                        }
                        doneBarrier.await();
                    }
                } catch (Exception error) {
                    throw new RuntimeException(error);
                }
            });
            workers[t].start();
        }
        try { warmupBarrier.await(); }
        catch (Exception error) { throw new RuntimeException(error); }
        final long start = System.nanoTime();
        for (int iteration = 0; iteration < iterations; ++iteration) {
            next.set(0);
            try { startBarrier.await(); doneBarrier.await(); }
            catch (Exception error) { throw new RuntimeException(error); }
        }
        for (Thread worker : workers) worker.join();
        final double seconds = (System.nanoTime() - start) * 1.0e-9;
        double checksum = 0.0; for (double value : values) checksum += value;
        System.out.printf("{\"mode\":\"gatk-java-logless\",\"pairs\":%d,\"threads\":%d,\"iterations\":%d,\"seconds\":%.12g,\"pairs_per_second\":%.12g,\"checksum\":%.12g}%n", inputs.size(), threads, iterations, seconds, inputs.size() * (double)iterations / seconds, checksum);
        if (!valuesOut.isEmpty()) try (BufferedWriter writer = new BufferedWriter(new FileWriter(valuesOut))) {
            for (double value : values) writer.write(Long.toHexString(Double.doubleToRawLongBits(value)) + "\n");
        }
    }

    private static void computeNativeBlocks(PairHMMNativeBinding binding, List<Pair> inputs, double[] values, int block) {
        for (int base = 0; base < inputs.size(); base += block) {
            final int n = Math.min(block, inputs.size() - base);
            final ReadDataHolder[] reads = new ReadDataHolder[n];
            final HaplotypeDataHolder[] haplotypes = new HaplotypeDataHolder[n];
            for (int i = 0; i < n; ++i) {
                Pair p = inputs.get(base + i);
                reads[i] = new ReadDataHolder(); reads[i].readBases = p.read; reads[i].readQuals = p.readQual;
                reads[i].insertionGOP = p.ins; reads[i].deletionGOP = p.del; reads[i].overallGCP = p.gcp;
                haplotypes[i] = new HaplotypeDataHolder(); haplotypes[i].haplotypeBases = p.hap;
            }
            final double[] matrix = new double[n * n];
            binding.computeLikelihoods(reads, haplotypes, matrix);
            for (int i = 0; i < n; ++i) values[base + i] = matrix[i * n + i];
        }
    }

    private static void writeTables(String path) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(path))) {
            for (int q = 0; q < 256; ++q) {
                final double error = q <= QualityUtils.MAX_QUAL
                        ? QualityUtils.qualToErrorProb((byte)q)
                        : Math.pow(10.0, -0.1 * q);
                writer.write("q " + Long.toHexString(Double.doubleToRawLongBits(error)) + "\n");
            }
            for (int ins = 0; ins < 256; ++ins) {
                for (int del = 0; del < 256; ++del) {
                    writer.write("m " + ins + " " + del + " " + Long.toHexString(Double.doubleToRawLongBits(PairHMMModel.matchToMatchProb((byte)ins, (byte)del))) + "\n");
                }
            }
        }
    }
}
