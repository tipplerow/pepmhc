
package pepmhc.engine;

/**
 * Categorizes the binding of peptides to MHC molecules.
 */
public enum BinderType {
    STRONG(50.0, 0.5),
    WEAK(500.0, 2.0),
    UNBOUND(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);

    private final double affinityThreshold;
    private final double percentileThreshold;

    private BinderType(double affinityThreshold, double percentileThreshold) {
        this.affinityThreshold = affinityThreshold;
        this.percentileThreshold = percentileThreshold;
    }

    /**
     * Returns the binder type for a given affinity and percentile
     * rank.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     *
     * @param percentile the percentile rank of the binding affinity
     * (relative to other peptides of the same length binding to the
     * same allele).
     *
     * @return the binder type for the given affinity and percentile
     * rank.
     */
    public static BinderType instance(double affinity, double percentile) {
        //
        // Percentile rank takes precedence over absolute affinity...
        //
        if (percentile <= STRONG.percentileThreshold)
            return STRONG;

        if (percentile <= WEAK.percentileThreshold)
            return WEAK;

        if (affinity <= STRONG.affinityThreshold)
            return STRONG;

        if (affinity <= WEAK.affinityThreshold)
            return WEAK;

        return UNBOUND;
    }

    /**
     * Returns the absolute affinity threshold for this binder type.
     *
     * @return the absolute affinity threshold for this binder type.
     */
    public double getAffinityThreshold() {
        return affinityThreshold;
    }

    /**
     * Returns the percentile rank threshold for this binder type.
     *
     * @return the percentile rank threshold for this binder type.
     */
    public double getPercentileThreshold() {
        return percentileThreshold;
    }
}
