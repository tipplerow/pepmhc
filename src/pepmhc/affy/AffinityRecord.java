
package pepmhc.affy;

import jene.peptide.Peptide;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the result of a peptide-MHC affinity measurement or
 * prediction.
 */
public final class AffinityRecord extends BindRecord {
    /**
     * Creates a new affinity record with an unset percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     */
    public AffinityRecord(Peptide peptide, double affinity) {
        super(peptide, affinity);
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
    public AffinityRecord(Peptide peptide, double affinity, double percentile) {
        super(peptide, affinity, percentile);
    }

    /**
     * Returns the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     *
     * @return the binding affinity expressed as an IC50 concentration
     * in nanomolar units.
     */
    @Override public double getAffinity() {
        return strength;
    }
}
