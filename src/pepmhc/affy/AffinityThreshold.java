
package pepmhc.affy;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import jam.app.JamProperties;
import jam.math.Percentile;

import jene.peptide.Peptide;

/**
 * Defines a threshold for peptide-MHC binding in terms of absolute
 * affinity, percentile rank, or both.
 */
public final class AffinityThreshold {
    private Affinity affinityThreshold;
    private Percentile percentileThreshold;

    private static AffinityThreshold global = null;

    private AffinityThreshold(Affinity affinityThreshold, Percentile percentileThreshold) {
        this.affinityThreshold = affinityThreshold;
        this.percentileThreshold = percentileThreshold;
        validate();
    }

    private void validate() {
        if (!isAffinityThresholdSet() && !isPercentileThresholdSet())
            throw new IllegalStateException("At least one threshold must be set.");
    }

    /**
     * Name of the system property that defines a global affinity
     * threshold (as an IC50 concentration in nanomolar units).
     */
    public static final String AFFINITY_THRESHOLD_PROPERTY = "pepmhc.affinityThreshold";

    /**
     * Name of the system property that defines a global percentile
     * rank affinity threshold.
     */
    public static final String PERCENTILE_THRESHOLD_PROPERTY = "pepmhc.percentileThreshold";

    /**
     * The standard affinity threshold that is usually applied.
     */
    public static final Affinity AFFINITY_THRESHOLD_STANDARD = Affinity.valueOf(500.0);

    /**
     * The standard percentile threshold that is usually applied.
     */
    public static final Percentile PERCENTILE_THRESHOLD_STANDARD = Percentile.valueOf(2.0);

    /**
     * The standard affinity and percentile thresholds that are
     * usually applied.
     */
    public static final AffinityThreshold STANDARD =
        AffinityThreshold.create(AFFINITY_THRESHOLD_STANDARD,
                                 PERCENTILE_THRESHOLD_STANDARD);

    /**
     * Creates a new affinity threshold with fixed attributes.
     *
     * @param affinity the affinity threshold (may be {@code null} to
     * apply only a percentile threshold).
     *
     * @param percentile the percentile threshold (may be {@code null} to
     * apply only an affinity threshold).
     *
     * @return a new affinity threshold with the specfied parameters.
     *
     * @throws RuntimeException if both thresholds are {@code null}.
     */
    public static AffinityThreshold create(Affinity affinity, Percentile percentile) {
        return new AffinityThreshold(affinity, percentile);
    }

    /**
     * Returns the global affinity threshold defined by system
     * properties.
     *
     * @return the global affinity threshold defined by system
     * properties.
     */
    public static AffinityThreshold global() {
        if (global == null)
            resolveGlobal();

        return global;
    }

    private static void resolveGlobal() {
        Affinity affinity = null;
        Percentile percentile = null;

        if (JamProperties.isSet(AFFINITY_THRESHOLD_PROPERTY))
            affinity = Affinity.valueOf(JamProperties.getRequiredDouble(AFFINITY_THRESHOLD_PROPERTY));

        if (JamProperties.isSet(PERCENTILE_THRESHOLD_PROPERTY))
            percentile = Percentile.valueOf(JamProperties.getRequiredDouble(PERCENTILE_THRESHOLD_PROPERTY));

        global = new AffinityThreshold(affinity, percentile);
    }

    /**
     * Counts the number of peptides that bind to an MHC molecule by
     * metrics of this threshold.
     *
     * @param records the peptide-MHC affinity records to examine.
     *
     * @return the number of peptides from the input collection that
     * are bound by the metrics of this threshold.
     */
    public int countBinders(Collection<AffinityRecord> records) {
        int result = 0;

        for (AffinityRecord record : records)
            if (isBound(record))
                ++result;

        return result;
    }

    /**
     * Finds the peptides that bind to an MHC molecule by metrics of
     * this threshold.
     *
     * @param records the peptide-MHC affinity records to examine.
     *
     * @return a set containing all peptides from the input collection
     * that are bound by the metrics of this threshold.
     */
    public Set<Peptide> getBinders(Collection<AffinityRecord> records) {
        Set<Peptide> binders = new HashSet<Peptide>();

        for (AffinityRecord record : records)
            if (isBound(record))
                binders.add(record.getPeptide());

        return binders;
    }


    /**
     * Determines whether a peptide is bound to an MHC molecule by the
     * metrics of this threshold.
     *
     * @param record a peptide-MHC affinity record to examine.
     *
     * @return {@code true} iff the peptide is bound by the metrics of
     * this threshold.
     */
    public boolean isBound(AffinityRecord record) {
        if (isAffinityThresholdSet()
            && record.getAffinity().LE(affinityThreshold))
            return true;

        if (isPercentileThresholdSet()
            && record.hasPercentile()
            && record.getPercentile().LE(percentileThreshold))
            return true;

        return false;
    }

    /**
     * Returns the absolute affinity threshold.
     *
     * @return the absolute affinity threshold.
     */
    public Affinity getAffinityThreshold() {
        return affinityThreshold;
    }

    /**
     * Returns the percentile rank threshold.
     *
     * @return the percentile rank threshold.
     */
    public Percentile getPercentileThreshold() {
        return percentileThreshold;
    }

    /**
     * Indicates whether the affinity threshold has been assigned.
     *
     * @return {@code true} iff the affinity threshold has been set.
     */
    public boolean isAffinityThresholdSet() {
        return affinityThreshold != null;
    }

    /**
     * Indicates whether the percentile rank threshold has been
     * assigned.
     *
     * @return {@code true} iff the percentile rank threshold has
     * been set.
     */
    public boolean isPercentileThresholdSet() {
        return percentileThreshold != null;
    }
}

