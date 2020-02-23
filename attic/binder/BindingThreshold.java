
package pepmhc.binder;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import jam.app.JamProperties;
import jam.peptide.Peptide;

/**
 * Defines a threshold for peptide-MHC binding in terms of absolute
 * affinity, percentile rank, or both.
 */
public final class BindingThreshold {
    private double affinityThreshold;
    private double percentileThreshold;

    private static BindingThreshold global = null;

    private BindingThreshold(double affinityThreshold, double percentileThreshold) {
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
    public static final String AFFINITY_THRESHOLD_PROPERTY = "pepmhc.binder.affinityThreshold";

    /**
     * Name of the system property that defines a global percentile
     * rank binding threshold.
     */
    public static final String PERCENTILE_THRESHOLD_PROPERTY = "pepmhc.binder.percentileThreshold";

    /**
     * Returns the global binding threshold defined by system
     * properties.
     *
     * @return the global binding threshold defined by system
     * properties.
     */
    public static BindingThreshold global() {
        if (global == null)
            resolveGlobal();

        return global;
    }

    private static void resolveGlobal() {
        double affinity = JamProperties.getOptionalDouble(AFFINITY_THRESHOLD_PROPERTY, Double.NaN);
        double percentile = JamProperties.getOptionalDouble(PERCENTILE_THRESHOLD_PROPERTY, Double.NaN);

        global = new BindingThreshold(affinity, percentile);
    }

    /**
     * Returns an object with a specified affinity threshold and an
     * unset percentile rank threshold.
     *
     * @param affinityThreshold the desired affinity threshold.
     *
     * @return an object with the specified affinity threshold and an
     * unset percentile rank threshold.
     */
    public static BindingThreshold forAffinity(double affinityThreshold) {
        return new BindingThreshold(affinityThreshold, Double.NaN);
    }

    /**
     * Returns an object with a specified percentile rank threshold
     * and an unset affinity threshold.
     *
     * @param percentileThreshold the desired percentile threshold.
     *
     * @return an object with the specified percentile rank threshold
     * and an unset affinity threshold.
     */
    public static BindingThreshold forPercentile(double percentileThreshold) {
        return new BindingThreshold(Double.NaN, percentileThreshold);
    }

    /**
     * Returns a new object with the percentile rank threshold of this
     * object and a specified affinity threshold.
     *
     * @param affinityThreshold the desired affinity threshold.
     *
     * @return a new object with the percentile rank threshold of this
     * object and the specified affinity threshold.
     */
    public BindingThreshold andAffinity(double affinityThreshold) {
        if (isAffinityThresholdSet())
            throw new IllegalStateException("Affinity threshold is already set.");

        return new BindingThreshold(affinityThreshold, this.percentileThreshold);
    }

    /**
     * Returns a new object with the affinity threshold of this object
     * and a specified percentile rank threshold.
     *
     * @param percentileThreshold the desired percentile threshold.
     *
     * @return a new object with the affinity threshold of this object
     * and the specified percentile rank threshold.
     */
    public BindingThreshold andPercentile(double percentileThreshold) {
        if (isPercentileThresholdSet())
            throw new IllegalStateException("Percentile threshold is already set.");

        return new BindingThreshold(this.affinityThreshold, percentileThreshold);
    }

    /**
     * Counts the number of peptides that bind to an MHC molecule by
     * metrics of this threshold.
     *
     * @param records the peptide-MHC binding records to examine.
     *
     * @return the number of peptides from the input collection that
     * are bound by the metrics of this threshold.
     */
    public int countBinders(Collection<BindingRecord> records) {
        int result = 0;

        for (BindingRecord record : records)
            if (isBound(record))
                ++result;

        return result;
    }

    /**
     * Finds the peptides that bind to an MHC molecule by metrics of
     * this threshold.
     *
     * @param records the peptide-MHC binding records to examine.
     *
     * @return a set containing all peptides from the input collection
     * that are bound by the metrics of this threshold.
     */
    public Set<Peptide> getBinders(Collection<BindingRecord> records) {
        Set<Peptide> binders = new HashSet<Peptide>();

        for (BindingRecord record : records)
            if (isBound(record))
                binders.add(record.getPeptide());

        return binders;
    }


    /**
     * Determines whether a peptide is bound to an MHC molecule by the
     * metrics of this threshold.
     *
     * @param record a peptide-MHC binding record to examine.
     *
     * @return {@code true} iff the peptide is bound by the metrics of
     * this threshold.
     */
    public boolean isBound(BindingRecord record) {
        if (isAffinityThresholdSet()
            && record.isAffinitySet()
            && record.getAffinity() <= affinityThreshold)
            return true;

        if (isPercentileThresholdSet()
            && record.isPercentileSet()
            && record.getPercentile() <= percentileThreshold)
            return true;

        return false;
    }

    /**
     * Returns the absolute affinity threshold.
     *
     * @return the absolute affinity threshold.
     */
    public double getAffinityThreshold() {
        return affinityThreshold;
    }

    /**
     * Returns the percentile rank threshold.
     *
     * @return the percentile rank threshold.
     */
    public double getPercentileThreshold() {
        return percentileThreshold;
    }

    /**
     * Indicates whether the affinity threshold has been assigned.
     *
     * @return {@code true} iff the affinity threshold has been set.
     */
    public boolean isAffinityThresholdSet() {
        return !Double.isNaN(affinityThreshold);
    }

    /**
     * Indicates whether the percentile rank threshold has been
     * assigned.
     *
     * @return {@code true} iff the percentile rank threshold has
     * been set.
     */
    public boolean isPercentileThresholdSet() {
        return !Double.isNaN(percentileThreshold);
    }
}

