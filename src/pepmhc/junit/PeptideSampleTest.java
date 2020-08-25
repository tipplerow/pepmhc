
package pepmhc.junit;

import java.util.Set;
import jean.peptide.Peptide;
import pepmhc.calc.PeptideSample;

import org.junit.*;
import static org.junit.Assert.*;

public class PeptideSampleTest {
    static {
        System.setProperty(PeptideSample.PEPTIDE_FILE_PROPERTY, "data/test/peptide_flat.txt");
        System.setProperty(PeptideSample.SAMPLE_COUNT_PROPERTY, "2");
        System.setProperty(PeptideSample.SAMPLE_SIZE_PROPERTY,  "4");
    }

    @Test public void testGlobal() {
        Set<Peptide> sample0 = PeptideSample.global().viewSample(0);
        Set<Peptide> sample1 = PeptideSample.global().viewSample(1);

        assertEquals(4, sample0.size());
        assertEquals(4, sample1.size());

        for (Peptide peptide : sample0)
            assertFalse(sample1.contains(peptide));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.PeptideSampleTest");
    }
}
