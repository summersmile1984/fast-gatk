import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.*;
import java.util.*;

public class P2 {
    static Allele al(String name, boolean ref) {
        if (name.equals("*")) return Allele.create("*", ref);
        return Allele.create(name, ref);
    }
    public static void main(String[] a) {
        // 1. what does getBaseString/toString give for a symbolic allele?
        Allele sym = Allele.create("<NON_REF>", false);
        System.out.println("symbolic: base=" + sym.getBaseString() + " display=" + sym.getDisplayString()
            + " toString=" + sym.toString() + " len=" + sym.length()
            + " baseslen=" + sym.getBases().length + " isSymbolic=" + sym.isSymbolic()
            + " equals(Allele.NON_REF_ALLELE)=" + sym.equals(Allele.NON_REF_ALLELE));
        System.out.println("SPAN_DEL: " + Allele.SPAN_DEL.getBaseString() + " len=" + Allele.SPAN_DEL.length()
            + " isSymbolic=" + Allele.SPAN_DEL.isSymbolic() + " isNonRef=" + Allele.SPAN_DEL.isNonReference());
        // 2. END attribute survival across the trim
        List<Allele> alleles = Arrays.asList(Allele.create("AAAA", true), Allele.create("AACA", false));
        Genotype g = GenotypeBuilder.create("S1", Arrays.asList(alleles.get(0), alleles.get(1)));
        VariantContext vc = new VariantContextBuilder("p", "chr1", 100, 103, alleles)
            .attribute("END", 103).attribute("DP", 20).genotypes(Collections.singletonList(g)).make();
        VariantContext out = org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils.reverseTrimAlleles(vc);
        System.out.println("END-case in : start=" + vc.getStart() + " end=" + vc.getEnd()
            + " ENDattr=" + vc.getAttribute("END") + " attrs=" + vc.getAttributes());
        System.out.println("END-case out: start=" + out.getStart() + " end=" + out.getEnd()
            + " ENDattr=" + out.getAttribute("END") + " attrs=" + out.getAttributes()
            + " alleles=" + out.getAlleles());
        // 3. encode the trimmed record through VCFEncoder with an END-bearing header
        VCFHeader header = new VCFHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine("END", 1, VCFHeaderLineType.Integer, "End"));
        header.addMetaDataLine(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "GT"));
        header.addMetaDataLine(new VCFContigHeaderLine(Collections.singletonMap("ID", "chr1"), 0));
        header.addMetaDataLine(new VCFInfoHeaderLine("AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "AC"));
        VariantContext withEnd = new VariantContextBuilder(out).attribute("END", 103).make();
        VCFEncoder enc = new VCFEncoder(header, false, false);
        System.out.println("encoded-with-END-attr : " + enc.encode(withEnd).trim());
        VariantContext noEnd = new VariantContextBuilder(out).rmAttribute("END").make();
        System.out.println("encoded-no-END-attr  : " + enc.encode(noEnd).trim());
    }
}
