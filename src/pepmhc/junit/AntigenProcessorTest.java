
package pepmhc.junit;

import java.util.Arrays;
import java.util.List;

import jam.junit.NumericTestBase;
import jam.peptide.Peptide;

import pepmhc.agpro.AntigenProcessor;
import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

import org.junit.*;
import static org.junit.Assert.*;

public class AntigenProcessorTest extends NumericTestBase {
    private static final Peptide A3GALT2 =
        Peptide.parse("MALKEGLRAWKRIFWRQILLTLGLLGLFLYGLPKFRHLEALIPMGVCPSATMSQLRDNFT" +
                      "GALRPWARPEVLTCTPWGAPIIWDGSFDPDVAKQEARQQNLTIGLTIFAVGRYLEKYLER" +
                      "FLETAEQHFMAGQSVMYYVFTELPGAVPRVALGPGRRLPVERVARERRWQDVSMARMRTL" +
                      "HAALGGLPGREAHFMFCMDVDQHFSGTFGPEALAESVAQLHSWHYHWPSWLLPFERDAHS" +
                      "AAAMAWGQGDFYNHAAVFGGSVAALRGLTAHCAGGLDWDRARGLEARWHDESHLNKFFWL" +
                      "HKPAKVLSPEFCWSPDIGPRAEIRRPRLLWAPKGYRLLRN");

    private static final Peptide BEX4 =
        Peptide.parse("MESKEELAANNLNGENAQQENEGGEQAPTQNEEESRHLGGGEGQKPGGNIRRGRVRRLVP" +
                      "NFRWAIPNRHIEHNEARDDVERFVGQMMEIKRKTREQQMRHYMRFQTPEPDNHYDFCLIP");

    @Test public void testNetChopTAPLenient() {
        AntigenProcessor processor =
            AntigenProcessor.resolve("data/test/NetChop_TAP_Lenient.prop");

        NetChop netChop = processor.getNetChop();

        assertTrue(Arrays.equals(new int[] { 9, 10 }, netChop.getLengths()));
        assertDouble(0.5, netChop.getThreshold());

        TAP tap = processor.getTAP();

        assertDouble(0.2, tap.getAlpha());
        assertDouble(1.0, tap.getThreshold());

        if (!NetChop.isInstalled())
            return;

        List<Peptide> fragments = processor.process(BEX4);

        assertEquals(List.of(Peptide.parse("HLGGGEGQK"),
                             Peptide.parse("AIPNRHIEH"),
                             Peptide.parse("KRKTREQQM"),
                             Peptide.parse("QTPEPDNHY"),
                             Peptide.parse("PNRHIEHNEA"),
                             Peptide.parse("FQTPEPDNHY")),
                     fragments);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AntigenProcessorTest");
    }
}
