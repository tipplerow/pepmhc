
package pepmhc.stab;

import java.util.Collection;
import java.util.HashMap;

import jam.util.IndexView;

import jean.peptide.Peptide;

/**
 * Provides a read-only view of stability records indexed by peptide.
 */
public final class StabilityMap extends IndexView<Peptide, StabilityRecord> {
    private StabilityMap(Collection<StabilityRecord> records) {
        super(new HashMap<Peptide, StabilityRecord>(), records, x -> x.getPeptide());
    }

    /**
     * Creates a new stability map.
     *
     * @param records the records to index.
     *
     * @return a new stability map containing the specified records.
     */
    public static StabilityMap create(Collection<StabilityRecord> records) {
        return new StabilityMap(records);
    }
}
