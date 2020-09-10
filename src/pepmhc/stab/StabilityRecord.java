
package pepmhc.stab;

import jam.math.Percentile;

import jene.chem.HalfLife;
import jene.peptide.Peptide;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the result of a peptide-MHC stability measurement or
 * prediction.
 */
public final class StabilityRecord extends BindRecord {
    private final HalfLife halfLife;

    /**
     * Creates a new stability record with an unset percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param halfLife the half-life of the MHC-bound state (hours).
     */
    public StabilityRecord(Peptide peptide, HalfLife halfLife) {
        this(peptide, halfLife, null);
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
    public StabilityRecord(Peptide peptide, HalfLife halfLife, Percentile percentile) {
        super(peptide, percentile);
        this.halfLife = halfLife;
    }

    /**
     * Returns the half-life of the MHC-bound state (in hours).
     *
     * @return the half-life of the MHC-bound state (in hours).
     */
    @Override public HalfLife getHalfLife() {
        return halfLife;
    }

    /**
     * Returns the binding strength (the half-life).
     *
     * @return the binding strength (the half-life).
     */
    @Override public double getStrength() {
        return halfLife.doubleValue();
    }
}
