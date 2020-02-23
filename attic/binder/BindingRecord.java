
package pepmhc.binder;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import jam.math.DoubleRange;
import jam.peptide.Peptide;

/**
 * Encapsulates the result of a peptide-MHC binding measurement or
 * prediction.
 */
public final class BindingRecord {
    private final Peptide peptide;
    private final double  affinity;
    private final double  percentile;

    /**
     * Creates a new binding record with an unset percentile rank.
     *
     * @param peptide the binding peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     */
    public BindingRecord(Peptide peptide, double affinity) {
        this(peptide, affinity, Double.NaN);
    }

    /**
     * Creates a new binding record.
     *
     * @param peptide the binding peptide.
     *
     * @param affinity the binding affinity expressed as an IC50
     * concentration in nanomolar units.
     *
     * @param percentile the percentile rank of the binding affinity
     * (relative to other peptides of the same length binding to the
     * same allele).
     */
    public BindingRecord(Peptide peptide, double affinity, double percentile) {
        this.peptide = peptide;
        this.affinity = affinity;
        this.percentile = percentile;

        validate();
    }

    private void validate() {
        DoubleRange.NON_NEGATIVE.validate("binding affinity", affinity);
        DoubleRange.PERCENTILE.validate("percentile rank", percentile);
    }

    public static Map<Peptide, BindingRecord> map(Collection<BindingRecord> records) {
        Map<Peptide, BindingRecord> map =
            new HashMap<Peptide, BindingRecord>(records.size());

        for (BindingRecord record : records)
            map.put(record.getPeptide(), record);

        return map;
    }

    public Peptide getPeptide() {
        return peptide;
    }

    public double getAffinity() {
        return affinity;
    }

    public double getPercentile() {
        return percentile;
    }

    public BinderType getBinderType() {
        return BinderType.instance(affinity, percentile);
    }

    public boolean isAffinitySet() {
        return !Double.isNaN(affinity);
    }

    public boolean isPercentileSet() {
        return !Double.isNaN(percentile);
    }

    @Override public String toString() {
        return String.format("BindingRecord(%s, %.1f, %.2f)", peptide.formatString(), affinity, percentile);
    }
}
