
package pepmhc.junit;

import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.stab.StabilityRecord;
import pepmhc.stab.StabilityPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class StabilityPredictorTest {
    private final Allele A0101 = Allele.instance("HLA-A*01:01");
    private final Allele A0201 = Allele.instance("HLA-A*02:01");

    @Test public void testExecutable() {
        System.out.print("NetStabPan ");
        System.out.print(StabilityPredictor.isInstalled() ? "IS" : "is NOT");
        System.out.print(" installed.");
        System.out.println();
    }

    @Test public void testPredict1() {
        if (!StabilityPredictor.isInstalled())
            return;

        assertEquals(0.43, StabilityPredictor.predict(A0101, Peptide.parse("AAAWYLWEV")).getHalfLife(), 0.01);
        assertEquals(0.11, StabilityPredictor.predict(A0101, Peptide.parse("AEFGPWQTV")).getHalfLife(), 0.01);

        assertEquals(7.61, StabilityPredictor.predict(A0201, Peptide.parse("AAAWYLWEV")).getHalfLife(), 0.01);
        assertEquals(0.50, StabilityPredictor.predict(A0201, Peptide.parse("AEFGPWQTV")).getHalfLife(), 0.01);

        Peptide pep8  = Peptide.parse("AAAWYLWE"); 
        Peptide pep9  = Peptide.parse("AAAWYLWEV");
        Peptide pep10 = Peptide.parse("AAAWYLWEVK");
        Peptide pep11 = Peptide.parse("AAAWYLWEVKL");

        List<Peptide> peptides = List.of(pep11, pep8, pep10, pep9);
        List<StabilityRecord> records = StabilityPredictor.predict(A0201, peptides);

        assertEquals(pep11, records.get(0).getPeptide());
        assertEquals(pep8,  records.get(1).getPeptide());
        assertEquals(pep10, records.get(2).getPeptide());
        assertEquals(pep9,  records.get(3).getPeptide());

        assertEquals(0.72, records.get(0).getHalfLife(), 0.01);
        assertEquals(0.13, records.get(1).getHalfLife(), 0.01);
        assertEquals(0.20, records.get(2).getHalfLife(), 0.01);
        assertEquals(7.61, records.get(3).getHalfLife(), 0.01);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.StabilityPredictorTest");
    }
}
