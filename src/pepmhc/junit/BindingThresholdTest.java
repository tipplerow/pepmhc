
package pepmhc.junit;

import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.binder.BindingThreshold;

import org.junit.*;
import static org.junit.Assert.*;

public class BindingThresholdTest {
    private static final Peptide peptide = Peptide.parse("AAA"); // Any peptide will do...

    private static final BindingRecord lowAffinity = new BindingRecord(peptide, 10.0,    99.9);
    private static final BindingRecord lowRank     = new BindingRecord(peptide, 10000.0,  0.1);
    private static final BindingRecord nonBinder   = new BindingRecord(peptide, 10000.0, 99.9);

    @Test public void testAndAffinity() {
        BindingThreshold threshold =
            BindingThreshold.forPercentile(2.0).andAffinity(500.0);

        assertTrue(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }
    
    @Test public void testAndPercentile() {
        BindingThreshold threshold =
            BindingThreshold.forAffinity(500.0).andPercentile(2.0);

        assertTrue(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }
    
    @Test public void testForAffinity() {
        BindingThreshold threshold = BindingThreshold.forAffinity(500.0);

        assertTrue(threshold.isAffinityThresholdSet());
        assertFalse(threshold.isPercentileThresholdSet());

        assertTrue(threshold.isBound(lowAffinity));
        assertFalse(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }

    @Test public void testForPercentile() {
        BindingThreshold threshold = BindingThreshold.forPercentile(2.0);

        assertFalse(threshold.isAffinityThresholdSet());
        assertTrue(threshold.isPercentileThresholdSet());

        assertFalse(threshold.isBound(lowAffinity));
        assertTrue(threshold.isBound(lowRank));
        assertFalse(threshold.isBound(nonBinder));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.BindingThresholdTest");
    }
}
