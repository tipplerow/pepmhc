
package pepmhc.junit;

import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.engine.Predictor;
import pepmhc.engine.net.NetStabPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class NetStabPredictorTest {
    private final Allele A0101 = Allele.instance("HLA-A*01:01");
    private final Allele A0201 = Allele.instance("HLA-A*02:01");
    private final Predictor predictor = NetStabPredictor.INSTANCE;

    @Test public void testExecutable() {
        System.out.print("NetStab ");
        System.out.print(predictor.isInstalled() ? "IS" : "is NOT");
        System.out.print(" installed.");
        System.out.println();
    }

    @Test public void testPredict1() {
        if (!predictor.isInstalled())
            return;

        assertEquals(0.43, predictor.predict(A0101, Peptide.parse("AAAWYLWEV")).getAffinity(), 0.01);
        assertEquals(0.11, predictor.predict(A0101, Peptide.parse("AEFGPWQTV")).getAffinity(), 0.01);

        assertEquals(7.61, predictor.predict(A0201, Peptide.parse("AAAWYLWEV")).getAffinity(), 0.01);
        assertEquals(0.50, predictor.predict(A0201, Peptide.parse("AEFGPWQTV")).getAffinity(), 0.01);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetStabPredictorTest");
    }
}
