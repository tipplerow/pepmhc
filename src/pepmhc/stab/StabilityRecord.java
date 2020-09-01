
package pepmhc.stab;

import jene.peptide.Peptide;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the result of a peptide-MHC stability measurement or
 * prediction.
 */
public final class StabilityRecord extends BindRecord {
    /**
     * Creates a new stability record with an unset percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param halfLife the half-life of the MHC-bound state (hours).
     */
    public StabilityRecord(Peptide peptide, double halfLife) {
        super(peptide, halfLife);
    }

    /**
     * Creates a new stability record.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param halfLife the half-life of the MHC-bound state (hours).
     *
     * @param percentile the percentile rank of the half-life
     * (relative to other peptides of the same length bound to
     * the same allele).
     */
    public StabilityRecord(Peptide peptide, double halfLife, double percentile) {
        super(peptide, halfLife, percentile);
    }

    /**
     * Returns the half-life of the MHC-bound state (in hours).
     *
     * @return the half-life of the MHC-bound state (in hours).
     */
    @Override public double getHalfLife() {
        return strength;
    }
}
