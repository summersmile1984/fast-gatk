import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Small release-pinned oracle for groupEvidence-before-marginalize order. */
public final class FragmentAggregationGatkOracle {
    private static final class Evidence implements Locatable {
        final String fragment;
        Evidence(final String fragment) { this.fragment = fragment; }
        public String getContig() { return "1"; }
        public int getStart() { return 1; }
        public int getEnd() { return 1; }
    }

    private static String bits(final double value) {
        return String.format("%016x", Double.doubleToRawLongBits(value));
    }

    public static void main(final String[] args) {
        final Allele ref0 = Allele.create("A", true);
        final Allele ref1 = Allele.create("C", false);
        final Allele alt = Allele.create("G", false);
        final List<Evidence> evidence = Arrays.asList(
                new Evidence("A"), new Evidence("B"),
                new Evidence("A"), new Evidence("B"));
        final double ninf = Double.NEGATIVE_INFINITY;
        final double[][][] values = new double[][][]{{
                {-1.0, -2.0, ninf, -3.0},
                {-2.0, -4.0, -3.0, -5.0},
                {-0.25, ninf, -0.75, -0.5}}};
        final AlleleLikelihoods<Evidence, Allele> reads =
                AlleleLikelihoods.createAlleleLikelihoods(
                        new IndexedAlleleList<>(ref0, ref1, alt),
                        new IndexedSampleList("S"), List.of(evidence), null, values);
        final AlleleLikelihoods<Evidence, Allele> fragments = reads.groupEvidence(
                item -> item.fragment,
                group -> new Evidence(group.get(0).fragment));
        final Map<Allele, List<Allele>> alleleMap = new LinkedHashMap<>();
        alleleMap.put(ref0, Arrays.asList(ref0, ref1));
        alleleMap.put(alt, List.of(alt));
        final AlleleLikelihoods<Evidence, Allele> alleles = fragments.marginalize(alleleMap);

        final LikelihoodMatrix<Evidence, Allele> fragmentMatrix = fragments.sampleMatrix(0);
        final String[] fragmentBits = new String[6];
        for (int evidenceIndex = 0; evidenceIndex < fragmentMatrix.evidenceCount(); ++evidenceIndex) {
            final int base = fragmentMatrix.getEvidence(evidenceIndex).fragment.equals("A") ? 0 : 3;
            for (int haplotype = 0; haplotype < 3; ++haplotype)
                fragmentBits[base + haplotype] = bits(fragmentMatrix.get(haplotype, evidenceIndex));
        }
        final LikelihoodMatrix<Evidence, Allele> alleleMatrix = alleles.sampleMatrix(0);
        final String[] alleleBits = new String[4];
        for (int evidenceIndex = 0; evidenceIndex < alleleMatrix.evidenceCount(); ++evidenceIndex) {
            final int base = alleleMatrix.getEvidence(evidenceIndex).fragment.equals("A") ? 0 : 2;
            alleleBits[base] = bits(alleleMatrix.get(0, evidenceIndex));
            alleleBits[base + 1] = bits(alleleMatrix.get(1, evidenceIndex));
        }
        System.out.println("{\"fragment_bits\":[\"" + String.join("\",\"", fragmentBits)
                + "\"],\"allele_bits\":[\"" + String.join("\",\"", alleleBits) + "\"]}");
    }
}
