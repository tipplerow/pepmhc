
package pepmhc.bind;

import jam.math.DoubleRange;

import jene.peptide.Peptide;

/**
 * Provides the base class for all binding records managed by the
 * {@code pepmhc} project.
 */
public abstract class BindRecord {
    /**
     * The peptide described by this record.
     */
    protected final Peptide peptide;

    /**
     * The binding strength (the affinity or stability metric).
     */
    protected final double strength;

    /**
     * The percentile rank of the binding strength (relative to other
     * peptides of the same length binding to the same allele).
     */
    protected final double percentile;

    /**
     * Creates a new record describing a given peptide with a fixed
     * binding strength (and an unset percentile rank).
     *
     * @param peptide the peptided described by the new record.
     *
     * @param strength the binding strength (the affinity or stability
     * metric).
     *
     * @throws RuntimeException unless the binding strength is positive.
     */
    protected BindRecord(Peptide peptide, double strength) {
        this(peptide, strength, Double.NaN);
    }

    /**
     * Creates a new record describing a given peptide with a fixed
     * binding strength and percentile rank.
     *
     * @param peptide the peptided described by the new record.
     *
     * @param strength the binding strength (the affinity or stability
     * metric).
     *
     * @param percentile the percentile rank of the binding strength.
     *
     * @throws RuntimeException unless the binding strength is positive
     * and the percentile is valid.
     */
    protected BindRecord(Peptide peptide, double strength, double percentile) {
        this.peptide = peptide;
        this.strength = strength;
        this.percentile = percentile;

        validate();
    }

    private void validate() {
        DoubleRange.NON_NEGATIVE.validate("binding strength", strength);

        if (!Double.isNaN(percentile))
            DoubleRange.PERCENTILE.validate("percentile rank", percentile);
    }

    /**
     * Returns the peptide described by this record.
     *
     * @return the peptide described by this record.
     */
    public final Peptide getPeptide() {
        return peptide;
    }

    /**
     * Returns the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     *
     * @return the binding affinity expressed as an IC50 concentration
     * in nanomolar units ({@code Double.NaN} unless this is an affinity
     * record).
     */
    public double getAffinity() {
        return Double.NaN;
    }

    /**
     * Returns the half-life of the MHC-bound state (in hours).
     *
     * @return the half-life of the MHC-bound state expressed in hours
     * ({@code Double.NaN} unless this is a stability record).
     */
    public double getHalfLife() {
        return Double.NaN;
    }

    /**
     * Returns the percentile rank of the binding strength (relative
     * to other peptides of the same length bound to the same allele).
     *
     * @return the percentile rank of the binding strength.
     */
    public final double getPercentile() {
        return percentile;
    }

    /**
     * Identifies records with a valid affinity.
     *
     * @return {@code false} iff the affinity is {@code Double.NaN}.
     */
    public boolean isAffinitySet() {
        return !Double.isNaN(getAffinity());
    }

    /**
     * Identifies records with a valid half-life.
     *
     * @return {@code false} iff the half-life is {@code Double.NaN}.
     */
    public boolean isHalfLifeSet() {
        return !Double.isNaN(getHalfLife());
    }

    /**
     * Identifies records with a valid percentile.
     *
     * @return {@code false} iff the percentile is {@code Double.NaN}.
     */
    public boolean isPercentileSet() {
        return !Double.isNaN(percentile);
    }

    @Override public String toString() {
        return String.format("%s(%s, %.2f, %.2f)", getClass().getSimpleName(),
                             peptide.formatString(), strength, percentile);
    }
}
