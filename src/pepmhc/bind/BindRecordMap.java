
package pepmhc.bind;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import jam.lang.JamException;
import jam.util.MapWrapper;

import jene.peptide.Peptide;

/**
 * Provides a map of binding records.
 */
public final class BindRecordMap extends MapWrapper<Peptide, BindRecord> {
    private BindRecordMap(Map<Peptide, BindRecord> map) {
        super(map);
    }

    /**
     * Creates a new binding record map backed by a {@code HashMap}.
     *
     * @return a new binding record map backed by a {@code HashMap}.
     */
    public static BindRecordMap hash() {
        return new BindRecordMap(new HashMap<Peptide, BindRecord>());
    }

    /**
     * Creates a new binding record map backed by a {@code HashMap}
     * and filled with a collection of binding records.
     *
     * @param records the records that will populate the new map.
     *
     * @return a new binding record map backed by a {@code HashMap} 
     * and filled with the specified binding records.
     */
    public static BindRecordMap hash(Collection<? extends BindRecord> records) {
        BindRecordMap map = hash();
        map.addAll(records);
        return map;
    }

    /**
     * Creates a new binding record map backed by a {@code TreeMap}.
     *
     * @return a new binding record map backed by a {@code TreeMap}.
     */
    public static BindRecordMap tree() {
        return new BindRecordMap(new TreeMap<Peptide, BindRecord>());
    }

    /**
     * Creates a new binding record map backed by a {@code TreeMap}
     * and filled with a collection of binding records.
     *
     * @param records the records that will populate the new map.
     *
     * @return a new binding record map backed by a {@code TreeMap} 
     * and filled with the specified binding records.
     */
    public static BindRecordMap tree(Collection<? extends BindRecord> records) {
        BindRecordMap map = tree();
        map.addAll(records);
        return map;
    }

    /**
     * Adds a binding record to this map.
     *
     * @param record the record to add.
     */
    public void add(BindRecord record) {
        map.put(record.getPeptide(), record);
    }

    /**
     * Adds binding records to this map.
     *
     * @param records the records to add.
     */
    public void addAll(Collection<? extends BindRecord> records) {
        for (BindRecord record : records)
            add(record);
    }

    /**
     * Identifies peptides contained in this map.
     *
     * @param peptide a peptide to examine.
     *
     * @return {@code true} iff this map contains the specified
     * peptide.
     */
    public boolean contains(Peptide peptide) {
        return map.containsKey(peptide);
    }

    /**
     * Retrieves a record from this map.
     *
     * @param peptide the peptide of the desired record.
     *
     * @return the record for the specified peptide.
     *
     * @throws RuntimeException unless this map contains a matching
     * record.
     */
    public BindRecord require(Peptide peptide) {
        BindRecord record = map.get(peptide);

        if (record != null)
            return record;
        else
            throw JamException.runtime("Missing binding record for [%s].", peptide);
    }

    @Override public String toString() {
        return map.values().toString();
    }
}
