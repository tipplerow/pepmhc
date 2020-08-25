
package pepmhc.junit;

import jean.peptide.Peptide;

import pepmhc.affy.AffinityRecord;
import pepmhc.affy.AffinityThreshold;

import org.junit.*;
import static org.junit.Assert.*;

public class AffinityThresholdTest {
    private static final Peptide peptide = Peptide.instance("AAA"); // Any peptide will do...

    private static final AffinityRecord lowAffinity = new AffinityRecord(peptide, 10.0,    99.9);
    private static final AffinityRecord lowRank     = new AffinityRecord(peptide, 10000.0,  0.1);
    private static final AffinityRecord nonBinder   = new AffinityRecord(peptide, 10000.0, 99.9);

    @Test public void testAndAffinity() {
        AffinityThreshold threshold =
            AffinityThreshold.forPercentile(2.0).andAffinity(500.0);

        assertTrue(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }
    
    @Test public void testAndPercentile() {
        AffinityThreshold threshold =
            AffinityThreshold.forAffinity(500.0).andPercentile(2.0);

        assertTrue(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }
    
    @Test public void testForAffinity() {
        AffinityThreshold threshold = AffinityThreshold.forAffinity(500.0);

        assertTrue(threshold.isAffinityThresholdSet());
        assertFalse(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertFalse(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }

    @Test public void testForPercentile() {
        AffinityThreshold threshold = AffinityThreshold.forPercentile(2.0);

        assertFalse(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertFalse(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AffinityThresholdTest");
    }
}
