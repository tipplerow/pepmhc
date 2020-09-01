
package pepmhc.junit;

import java.util.List;

import jam.math.JamRandom;

import jene.hla.Allele;
import jene.hla.Genotype;
import jene.peptide.Peptide;

import pepmhc.affy.AffinityCache;
import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityPredictor;
import pepmhc.affy.AffinityThreshold;
import pepmhc.affy.net.NetMHCPredictor;
import pepmhc.calc.PresentationRateCalculator;

import org.junit.*;
import static org.junit.Assert.*;

public class PresentationRateCalculatorTest {
    static {
        System.setProperty(JamRandom.SEED_PROPERTY, "123456789");
        System.setProperty(AffinityCache.CACHE_DIRECTORY_PROPERTY, "data/cache");
        System.setProperty(NetMHCPredictor.EXECUTABLE_PATH_PROPERTY, "/Users/scott/local/netMHC-4.0/netMHC");
    }

    private final Allele A0201 = Allele.instance("HLA-A*02:01");
    private final Allele B4002 = Allele.instance("HLA-B*40:02");
    private final AffinityMethod method = AffinityMethod.NET_MHC;
    private final List<Peptide> peptides = Peptide.newNative(9, 1000);
    private final AffinityThreshold threshold = AffinityThreshold.forAffinity(500.0);

    @Test public void testNetMHC() {
        if (!method.getPredictor().isInstalled())
            return;

        PresentationRateCalculator calculator =
            PresentationRateCalculator.instance(method, threshold, peptides);

        assertEquals(0.035, calculator.compute(A0201), 0.001);
        assertEquals(0.015, calculator.compute(B4002), 0.001);
        assertEquals(0.050, calculator.compute(Genotype.instance(A0201, B4002)), 0.001);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.PresentationRateCalculatorTest");
    }
}
