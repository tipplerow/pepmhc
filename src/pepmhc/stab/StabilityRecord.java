
package pepmhc.stab;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import jam.math.DoubleRange;
import jam.peptide.Peptide;
import jam.util.MapUtil;

/**
 * Encapsulates the result of a peptide-MHC stability measurement or
 * prediction.
 */
public final class StabilityRecord {
    private final Peptide peptide;
    private final double  halfLife;
    private final double  percentile;

    /**
     * Creates a new stability record with an unset percentile rank.
     *
     * @param peptide the MHC-bound peptide.
     *
     * @param halfLife the half-life of the MHC-bound state (hours).
     */
    public StabilityRecord(Peptide peptide, double halfLife) {
        this(peptide, halfLife, Double.NaN);
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
        this.peptide = peptide;
        this.halfLife = halfLife;
        this.percentile = percentile;

        validate();
    }

    private void validate() {
        DoubleRange.NON_NEGATIVE.validate("half-life", halfLife);
        DoubleRange.PERCENTILE.validate("percentile rank", percentile);
    }

    /**
     * Creates a map of stability records indexed by peptide.
     *
     * @param records the records to map.
     *
     * @return a map of stability records indexed by peptide.
     *
     * @throws RuntimeException if the collection contains a duplicate
     * peptide.
     */
    public static Map<Peptide, StabilityRecord> map(Collection<StabilityRecord> records) {
        Map<Peptide, StabilityRecord> map = new HashMap<Peptide, StabilityRecord>(records.size());

        for (StabilityRecord record : records)
            MapUtil.putUnique(map, record.getPeptide(), record);

        return map;
    }

    /**
     * Returns the MHC-bound peptide.
     *
     * @return the MHC-bound peptide.
     */
    public Peptide getPeptide() {
        return peptide;
    }

    /**
     * Returns the half-life of the MHC-bound state (hours).
     *
     * @return the half-life of the MHC-bound state (hours).
     */
    public double getHalfLife() {
        return halfLife;
    }

    /**
     * Returns the percentile rank of the half-life (relative to other
     * peptides of the same length bound to the same allele).
     *
     * @return the percentile rank of the half-life (relative to other
     * peptides of the same length bound to the same allele).
     */
    public double getPercentile() {
        return percentile;
    }

    @Override public String toString() {
        return String.format("StabilityRecord(%s, %.2f, %.2f)", peptide.formatString(), halfLife, percentile);
    }
}
