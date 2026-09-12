// Unit-level oracle: what does GATK 4.6.2.0's own
// GATKVariantContextUtils.reverseTrimAlleles() do to a hand-built
// VariantContext?  Compiled against the pinned GATK fat jar (which bundles
// htsjdk), so the answer is the shipped bytecode's behaviour, not a reading of
// the Java sources.
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

public class TrimProbe {
    static Allele allele(String name, boolean ref) {
        if (name.equals("*")) return ref ? Allele.create("*", true) : Allele.SPAN_DEL;
        if (name.startsWith("<")) return Allele.create(name, ref);
        return Allele.create(name, ref);
    }

    static void probe(String label, int start, String... names) {
        List<Allele> alleles = new ArrayList<>();
        for (int i = 0; i < names.length; i++) alleles.add(allele(names[i], i == 0));
        List<Genotype> gts = new ArrayList<>();
        gts.add(GenotypeBuilder.create("S1", Arrays.asList(alleles.get(0), alleles.get(0))));
        VariantContext vc = new VariantContextBuilder("probe", "chr1", start,
                start + alleles.get(0).length() - 1, alleles).genotypes(gts).make();
        VariantContext out = GATKVariantContextUtils.reverseTrimAlleles(vc);
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%-34s in=%-28s", label, String.join("/", names)));
        sb.append(String.format(" type=%-12s isVariant=%-5s", vc.getType(), vc.isVariant()));
        sb.append(" out=").append(String.join("/", out.getAlleles().stream().map(Allele::getBaseString).toList()));
        sb.append(String.format(" start=%d end=%d unchanged=%s", out.getStart(), out.getEnd(), out == vc));
        System.out.println(sb);
    }

    public static void main(String[] args) {
        System.out.println("# GATKVariantContextUtils.reverseTrimAlleles on synthetic records");
        probe("multibase-ref-nonref", 100, "AAAA", "<NON_REF>");
        probe("multibase-ref-nonref-8", 100, "ACGTACGT", "<NON_REF>");
        probe("multibase-ref-star", 100, "AA", "*");
        probe("multibase-ref-star-nonref", 100, "AA", "*", "<NON_REF>");
        probe("ref-alt-shared3", 100, "AAAA", "AA");
        probe("ref-alt-shared1", 100, "AAAA", "AACA");
        probe("ref-alt-noshared", 100, "AAAA", "AACC");
        probe("guard-one-base-alt", 100, "AAA", "A");
        probe("guard-star-only-alt", 100, "AAA", "*");
        probe("empty-guard-acgtacgt", 100, "ACGTACGT", "ACGT");
        probe("ref-lt-alt-insertion", 100, "AA", "AAAA");
        probe("ref-2alt-mixed-len", 100, "AAAAA", "AAA", "CAAAA");
        probe("ref-star-long", 100, "AAAAA", "*");
        probe("symbolic-only", 100, "ACGT", "<DEL>");
        probe("symbolic-plus-nonref", 100, "ACGT", "<DEL>", "<NON_REF>");
        probe("ref-alt-symbolic-shared", 100, "ACGTAC", "ACGT", "<NON_REF>");
        probe("star-ref-multibase-3", 100, "AAA", "*", "<NON_REF>");
        System.out.println("# htsjdk=" + VariantContext.class.getPackage().getImplementationVersion());
    }
}
