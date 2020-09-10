
package pepmhc.junit;

import jam.math.Percentile;

import jene.peptide.Peptide;

import pepmhc.affy.Affinity;
import pepmhc.affy.AffinityRecord;
import pepmhc.affy.AffinityThreshold;

import org.junit.*;
import static org.junit.Assert.*;

public class AffinityThresholdTest {
    private static final Peptide peptide = Peptide.instance("AAA"); // Any peptide will do...

    private static final AffinityRecord lowAffinity =
        new AffinityRecord(peptide, Affinity.valueOf(10.0), Percentile.valueOf(99.9));

    private static final AffinityRecord lowRank =
        new AffinityRecord(peptide, Affinity.valueOf(10000.0), Percentile.valueOf(0.1));

    private static final AffinityRecord nonBinder =
        new AffinityRecord(peptide, Affinity.valueOf(10000.0), Percentile.valueOf(99.9));

    @Test public void testStandard() {
        AffinityThreshold threshold =
            AffinityThreshold.STANDARD;

        assertTrue(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }
    
    @Test public void testAffinityOnly() {
        AffinityThreshold threshold = AffinityThreshold.create(Affinity.valueOf(500.0), null);

        assertTrue(threshold.isAffinityThresholdSet());
        assertFalse(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertFalse(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }

    @Test public void testPercentileOnly() {
        AffinityThreshold threshold = AffinityThreshold.create(null, Percentile.valueOf(2.0));

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
