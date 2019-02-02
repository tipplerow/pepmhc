
package pepmhc.junit;

import java.util.List;

import jam.hla.Allele;
import jam.hla.Genotype;
import jam.math.JamRandom;
import jam.peptide.Peptide;

import pepmhc.binder.BindingThreshold;
import pepmhc.calc.PresentationRateCalculator;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;
import pepmhc.engine.net.NetMHCPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class PresentationRateCalculatorTest {
    static {
        System.setProperty(JamRandom.SEED_PROPERTY, "123456789");
        System.setProperty(AffinityCache.CACHE_DIRECTORY_PROPERTY, "data/test");
        System.setProperty(NetMHCPredictor.EXECUTABLE_PATH_PROPERTY, "/Users/scott/local/netMHC-4.0/netMHC");
    }

    private final Allele A0201 = Allele.instance("HLA-A*02:01");
    private final Allele B4002 = Allele.instance("HLA-B*40:02");
    private final PredictionMethod method = PredictionMethod.NET_MHC;
    private final List<Peptide> peptides = Peptide.newNative(9, 1000);
    private final BindingThreshold threshold = BindingThreshold.forAffinity(500.0);

    @Test public void testNetMHC() {
        if (!Predictor.isInstalled(method))
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
