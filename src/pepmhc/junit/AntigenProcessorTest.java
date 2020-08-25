
package pepmhc.junit;

import java.util.Arrays;
import java.util.List;

import jam.junit.NumericTestBase;

import jean.peptide.Peptide;

import pepmhc.agpro.AntigenProcessor;
import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

import org.junit.*;
import static org.junit.Assert.*;

public class AntigenProcessorTest extends NumericTestBase {
    private static final Peptide A3GALT2 =
        Peptide.instance("MALKEGLRAWKRIFWRQILLTLGLLGLFLYGLPKFRHLEALIPMGVCPSATMSQLRDNFT" +
                         "GALRPWARPEVLTCTPWGAPIIWDGSFDPDVAKQEARQQNLTIGLTIFAVGRYLEKYLER" +
                         "FLETAEQHFMAGQSVMYYVFTELPGAVPRVALGPGRRLPVERVARERRWQDVSMARMRTL" +
                         "HAALGGLPGREAHFMFCMDVDQHFSGTFGPEALAESVAQLHSWHYHWPSWLLPFERDAHS" +
                         "AAAMAWGQGDFYNHAAVFGGSVAALRGLTAHCAGGLDWDRARGLEARWHDESHLNKFFWL" +
                         "HKPAKVLSPEFCWSPDIGPRAEIRRPRLLWAPKGYRLLRN");

    private static final Peptide BEX4 =
        Peptide.instance("MESKEELAANNLNGENAQQENEGGEQAPTQNEEESRHLGGGEGQKPGGNIRRGRVRRLVP" +
                         "NFRWAIPNRHIEHNEARDDVERFVGQMMEIKRKTREQQMRHYMRFQTPEPDNHYDFCLIP");

    @Test public void testNetChopTAPLenient() {
        AntigenProcessor processor =
            AntigenProcessor.resolve("data/test/NetChop_TAP_Lenient.prop");

        TAP tap = processor.getTAP();

        assertDouble(0.2, tap.getAlpha());
        assertDouble(1.0, tap.getThreshold());

        if (!NetChop.isInstalled())
            return;

        List<Peptide> fragments = processor.process(BEX4);

        assertEquals(List.of(Peptide.instance("HLGGGEGQK"),
                             Peptide.instance("AIPNRHIEH"),
                             Peptide.instance("KRKTREQQM"),
                             Peptide.instance("QTPEPDNHY"),
                             Peptide.instance("PNRHIEHNEA"),
                             Peptide.instance("FQTPEPDNHY")),
                     fragments);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AntigenProcessorTest");
    }
}
