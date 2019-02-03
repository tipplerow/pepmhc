
package pepmhc.junit;

import java.util.Set;
import jam.peptide.Peptide;
import pepmhc.cache.SampleSet;

import org.junit.*;
import static org.junit.Assert.*;

public class SampleSetTest {
    static {
        System.setProperty(SampleSet.PEPTIDE_FILE_PROPERTY, "data/test/peptide_flat.txt");
        System.setProperty(SampleSet.SAMPLE_COUNT_PROPERTY, "2");
        System.setProperty(SampleSet.SAMPLE_SIZE_PROPERTY,  "4");
        System.setProperty(SampleSet.WITH_REPLACEMENT_PROPERTY, "false");
    }

    @Test public void testAll() {
        Set<Peptide> sample0 = SampleSet.global().viewSample(0);
        Set<Peptide> sample1 = SampleSet.global().viewSample(1);

        assertEquals(4, sample0.size());
        assertEquals(4, sample1.size());

        for (Peptide peptide : sample0)
            assertFalse(sample1.contains(peptide));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.SampleSetTest");
    }
}
