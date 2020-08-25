
package pepmhc.junit;

import jam.junit.NumericTestBase;

import jean.hla.Allele;
import jean.peptide.Peptide;

import pepmhc.stab.StabilityRecord;
import pepmhc.stab.proxy.AffinityProxyModel;
import pepmhc.stab.proxy.AffinityProxyStore;

import org.junit.*;
import static org.junit.Assert.*;


public class AffinityProxyModelTest {
    private final Allele A0101 = Allele.instance("HLA-A*01:01");
    private final Allele A0201 = Allele.instance("HLA-A*02:01");

    @Test public void testInstance() {
        AffinityProxyModel model1 = AffinityProxyStore.get(A0101);
        AffinityProxyModel model2 = AffinityProxyStore.get(A0201);

        assertEquals( 2.614877, model1.getIntercept(),   0.000001);
        assertEquals(-0.436525, model1.getCoefficient(), 0.000001);
        assertEquals( 2.691525, model2.getIntercept(),   0.000001);
        assertEquals(-0.448037, model2.getCoefficient(), 0.000001);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AffinityProxyModelTest");
    }
}
