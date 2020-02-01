
package pepmhc.junit;

import jam.junit.NumericTestBase;
import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.stab.AffinityProxyModel;
import pepmhc.stab.StabilityRecord;

import org.junit.*;
import static org.junit.Assert.*;


public class AffinityProxyModelTest {
    private final Allele A0101 = Allele.instance("HLA-A*01:01");
    private final Allele A0201 = Allele.instance("HLA-A*02:01");

    @Test public void testInstance() {
        AffinityProxyModel model1 = AffinityProxyModel.instance(A0101);
        AffinityProxyModel model2 = AffinityProxyModel.instance(A0201);

        assertEquals( 2.593924, model1.getIntercept(),   0.000001);
        assertEquals(-0.434550, model1.getCoefficient(), 0.000001);
        assertEquals( 2.675812, model2.getIntercept(),   0.000001);
        assertEquals(-0.446402, model2.getCoefficient(), 0.000001);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AffinityProxyModelTest");
    }
}
