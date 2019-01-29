
package pepmhc.junit;

import jam.peptide.Peptide;

import pepmhc.engine.Predictor;
import pepmhc.engine.net.NetMHCPanPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class NetMHCPanPredictorTest {
    static {
        System.setProperty(NetMHCPanPredictor.EXECUTABLE_PATH_PROPERTY, "/Users/scott/local/netMHCpan-4.0/netMHCpan");
    }

    private final Predictor predictor = NetMHCPanPredictor.INSTANCE;

    @Test public void testExecutable() {
        System.out.print("NetMHCPan ");
        System.out.print(NetMHCPanPredictor.isInstalled() ? "IS" : "is NOT");
        System.out.print(" installed.");
        System.out.println();
    }

    @Test public void testPredict1() {
        if (!NetMHCPanPredictor.isInstalled())
            return;

        assertEquals( 8706.7, predictor.predict("HLA-A*02:01", Peptide.parse("AEFGPWQTV")).getAffinity(), 0.1);
        assertEquals(33803.5, predictor.predict("HLA-A*01:01", Peptide.parse("AEFGPWQTV")).getAffinity(), 0.1);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetMHCPanPredictorTest");
    }
}
