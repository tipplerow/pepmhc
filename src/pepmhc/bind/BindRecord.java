
package pepmhc.bind;

import jam.math.DoubleRange;
import jam.math.Percentile;

import jene.chem.HalfLife;
import jene.peptide.Peptide;

import pepmhc.affy.Affinity;

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
     * The percentile rank of the binding strength (relative to other
     * peptides of the same length binding to the same allele).
     */
    protected final Percentile percentile;

    /**
     * Creates a new record describing a given peptide with a missing
     * percentile rank.
     *
     * @param peptide the peptided described by the new record.
     */
    protected BindRecord(Peptide peptide) {
        this(peptide, null);
    }

    /**
     * Creates a new record describing a given peptide with a fixed
     * percentile binding rank.
     *
     * @param peptide the peptided described by the new record.
     *
     * @param percentile the percentile rank of the binding strength.
     */
    protected BindRecord(Peptide peptide, Percentile percentile) {
        this.peptide = peptide;
        this.percentile = percentile;
    }

    /**
     * Returns the binding strength (the affinity or stability metric).
     *
     * @return the binding strength (the affinity or stability metric).
     */
    public abstract double getStrength();

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
     * in nanomolar units.
     *
     * @throws UnsupportedOperationException unless this is an affinity
     * record.
     */
    public Affinity getAffinity() {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns the half-life of the MHC-bound state (in hours).
     *
     * @return the half-life of the MHC-bound state expressed in hours.
     *
     * @throws UnsupportedOperationException unless this is a stability
     * record.
     */
    public HalfLife getHalfLife() {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns the percentile rank of the binding strength (relative
     * to other peptides of the same length bound to the same allele).
     *
     * @return the percentile rank of the binding strength ({@code null}
     * if the rank is unknown).
     */
    public final Percentile getPercentile() {
        return percentile;
    }

    /**
     * Identifies records with a valid percentile.
     *
     * @return {@code true} iff this record has a non-{@code null}
     * percentile.
     */
    public boolean hasPercentile() {
        return percentile != null;
    }

    @Override public String toString() {
        return String.format("%s(%s, %.2f, %.2f)", getClass().getSimpleName(),
                             peptide.formatString(), getStrength(), percentile.doubleValue());
    }
}
