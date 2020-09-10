
package pepmhc.affy;

import jam.math.Percentile;

import jene.peptide.Peptide;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the result of a peptide-MHC affinity measurement or
 * prediction.
 */
public final class AffinityRecord extends BindRecord {
    private final Affinity affinity;

    /**
     * Creates a new affinity record with a missing percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     */
    public AffinityRecord(Peptide peptide, Affinity affinity) {
        this(peptide, affinity, null);
    }

    /**
     * Creates a new affinity record.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     *
     * @param percentile the percentile rank of the binding affinity
     * (relative to other peptides of the same length binding to the
     * same allele).
     */
    public AffinityRecord(Peptide peptide, Affinity affinity, Percentile percentile) {
        super(peptide, percentile);
        this.affinity = affinity;
    }

    /**
     * Returns the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     *
     * @return the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     */
    @Override public Affinity getAffinity() {
        return affinity;
    }

    /**
     * Returns the binding strength (the affinity).
     *
     * @return the binding strength (the affinity).
     */
    @Override public double getStrength() {
        return affinity.doubleValue();
    }
}
