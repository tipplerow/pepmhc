
package pepmhc.stab;

import java.util.ArrayList;
import java.util.List;

import jam.math.DoubleRange;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;

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
     * Converts a binding record into a stability record.
     *
     * @param bindingRecord a binding record that actually has a
     * half-life in its affinity attribute.
     *
     * @return the stability record corresponding to the given binding
     * record.
     */
    public static StabilityRecord convert(BindingRecord bindingRecord) {
        return new StabilityRecord(bindingRecord.getPeptide(),
                                   bindingRecord.getAffinity(),
                                   bindingRecord.getPercentile());
    }

    /**
     * Converts a list of binding records into stability records.
     *
     * @param bindingRecords a list of binding records that actually
     * have half-lives in their affinity attributes.
     *
     * @return the stability records corresponding to the given
     * binding records.
     */
    public static List<StabilityRecord> convert(List<BindingRecord> bindingRecords) {
        List<StabilityRecord> stabilityRecords = new ArrayList<StabilityRecord>(bindingRecords.size());

        for (BindingRecord bindingRecord : bindingRecords)
            stabilityRecords.add(convert(bindingRecord));

        return stabilityRecords;
    }

    public Peptide getPeptide() {
        return peptide;
    }

    public double getHalfLife() {
        return halfLife;
    }

    public double getPercentile() {
        return percentile;
    }

    public boolean isHalfLifeSet() {
        return !Double.isNaN(halfLife);
    }

    public boolean isPercentileSet() {
        return !Double.isNaN(percentile);
    }

    @Override public String toString() {
        return String.format("StabilityRecord(%s, %.2f, %.2f)", peptide.formatString(), halfLife, percentile);
    }
}
