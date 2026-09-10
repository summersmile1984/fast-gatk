package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.read.FlowBasedKeyCodec;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.util.Arrays;

/**
 * Emits a small deterministic FlowBasedPairHMM corpus for the native Kokkos
 * recurrence.  It intentionally calls the same protected GATK method used by
 * FlowBasedHMMEngine, rather than reimplementing the recurrence in the test.
 */
public final class FlowPairHmmGatkOracle {
    private static final class ExposedPairHMM extends FlowBasedPairHMM {
        double score(final int start, final int[] hapKey, final byte[] hapOrder,
                     final int[] readKey, final byte[] readOrder,
                     final FlowBasedRead read, final byte[] insertion,
                     final byte[] deletion, final byte[] continuation) {
            initialize(readKey.length, hapKey.length);
            return subComputekeayLikelihoodGivenHaplotypeKeysLog10(
                    start, hapKey, hapOrder, readKey, readOrder, read,
                    insertion, deletion, continuation, true);
        }
    }

    private static FlowBasedRead read(final String sequence, final String flowOrder) {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("fastgatk-flow-oracle");
        record.setReadString(sequence);
        final byte[] qualities = new byte[sequence.length()];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setCigarString(sequence.length() + "M");
        final byte[] tp = new byte[sequence.length()];
        record.setAttribute("tp", tp);
        final FlowBasedArgumentCollection arguments = new FlowBasedArgumentCollection();
        arguments.keepBoundaryFlows = true;
        arguments.fillingValue = 0.001;
        final FlowBasedRead result = new FlowBasedRead(record, flowOrder, 12, arguments);
        result.applyAlignment();
        return result;
    }

    private static int firstFlow(final byte[] hapOrder, final byte readFirst) {
        for (int index = 0; index < hapOrder.length; ++index)
            if (hapOrder[index] == readFirst) return index;
        return 0;
    }

    private static void emit(final String name, final String readBases,
                             final String haplotypeBases, final String flowOrder) {
        final FlowBasedRead read = read(readBases, flowOrder);
        final byte[] hapBases = haplotypeBases.getBytes();
        final int[] hapKey = FlowBasedKeyCodec.baseArrayToKey(hapBases, flowOrder);
        final byte[] hapOrder = FlowBasedKeyCodec.getFlowToBase(flowOrder, hapKey.length);
        final int[] readKey = read.getKey();
        final byte[] readOrder = read.getFlowOrderArray();
        final byte[] insertion = new byte[readKey.length];
        final byte[] deletion = new byte[readKey.length];
        final byte[] continuation = new byte[readKey.length];
        Arrays.fill(insertion, (byte) 40);
        Arrays.fill(deletion, (byte) 40);
        Arrays.fill(continuation, (byte) 10);
        final double score = new ExposedPairHMM().score(
                firstFlow(hapOrder, readOrder[0]), hapKey, hapOrder,
                readKey, readOrder, read, insertion, deletion, continuation);
        System.out.println(name + "\t" + score + "\t0x" +
                Long.toHexString(Double.doubleToRawLongBits(score)) + "\t" +
                Arrays.toString(readKey) + "\t" + Arrays.toString(hapKey));
    }

    public static void main(final String[] args) {
        emit("match", "ACGTACGT", "ACGTACGT", "ACGT");
        emit("mismatch", "ACGTACGT", "AGGTACGT", "ACGT");
        emit("homopolymer", "AAACCCGGGTTT", "AACCCGGTTT", "ACGT");
        emit("phase", "TACGTACG", "TACGTTACG", "TGCA");
    }
}
