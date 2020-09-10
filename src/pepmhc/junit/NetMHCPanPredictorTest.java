
package pepmhc.junit;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.affy.AffinityPredictor;
import pepmhc.affy.net.NetMHCPanPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class NetMHCPanPredictorTest {
    private final Allele A0101 = Allele.instance("HLA-A*01:01");
    private final Allele A0201 = Allele.instance("HLA-A*02:01");
    private final AffinityPredictor predictor = NetMHCPanPredictor.INSTANCE;

    @Test public void testExecutable() {
        System.out.print("NetMHCPan ");
        System.out.print(predictor.isInstalled() ? "IS" : "is NOT");
        System.out.print(" installed.");
        System.out.println();
    }

    @Test public void testPredict1() {
        if (!predictor.isInstalled())
            return;

        assertEquals( 8706.7, predictor.predict(A0201, Peptide.instance("AEFGPWQTV")).getStrength(), 0.1);
        assertEquals(33803.5, predictor.predict(A0101, Peptide.instance("AEFGPWQTV")).getAffinity().doubleValue(), 0.1);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetMHCPanPredictorTest");
    }
}
