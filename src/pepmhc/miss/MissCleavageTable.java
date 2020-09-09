
package pepmhc.miss;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import jam.app.JamLogger;
import jam.io.LineReader;
import jam.util.PairKeyTable;

import jene.hugo.HugoSymbol;
import jene.tcga.TumorBarcode;

/**
 * Indexes missense-chop records by tumor barcode and HUGO symbol.
 */
public final class MissCleavageTable {
    //
    // Class-specific container...
    //
    private static final class RecordList extends ArrayList<MissCleavageRecord> {}

    // All records indexed by barcode (outer) and symbol (inner)...
    private final PairKeyTable<TumorBarcode, HugoSymbol, RecordList> table = PairKeyTable.hash();

    // Total number of records in the table...
    private int count = 0;

    private MissCleavageTable(Collection<MissCleavageRecord> records) {
        fillMap(records);
    }

    private void fillMap(Collection<MissCleavageRecord> records) {
        for (MissCleavageRecord record : records)
            addRecord(record);
    }

    private void addRecord(MissCleavageRecord record) {
        HugoSymbol   symbol  = record.getHugoSymbol();
        TumorBarcode barcode = record.getTumorBarcode();

        ++count;
        recordList(barcode, symbol).add(record);
    }

    private RecordList recordList(TumorBarcode barcode, HugoSymbol symbol) {
        RecordList recordList = table.get(barcode, symbol);

        if (recordList == null) {
            recordList = new RecordList();
            table.put(barcode, symbol, recordList);
        }

        return recordList;
    }

    /**
     * Populates a table by reading all records from a given file.
     *
     * @param fileName the path to the missense mutation file.
     *
     * @return a table containing all records in the given file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains properly formatted records.
     */
    public static MissCleavageTable load(String fileName) {
        return load(new File(fileName));
    }

    /**
     * Populates a table by reading all records from a given file.
     *
     * @param file the path to the missense mutation file.
     *
     * @return a table containing all records in the given file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains properly formatted records.
     */
    public static MissCleavageTable load(File file) {
        List<MissCleavageRecord> records =
            new ArrayList<MissCleavageRecord>();

        try (LineReader reader = LineReader.open(file)) {
            // Skip header line...
            reader.next();

            for (String line : reader)
                records.add(MissCleavageRecord.parse(line));
        }

        JamLogger.info("MissCleavageTable: Loaded [%d] records.", records.size());
        return load(records);
    }

    /**
     * Populates a table from a collection of records.
     *
     * @param records the records to be indexed in the table.
     *
     * @return a table containing all records in the given collection.
     */
    public static MissCleavageTable load(Collection<MissCleavageRecord> records) {
        return new MissCleavageTable(records);
    }

    /**
     * Identifies tumor-gene pairs contained in this mutation table.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @param symbol the gene of interest.
     *
     * @return {@code true} iff this table contains records for the
     * specified tumor-gene pair.
     */
    public boolean contains(TumorBarcode barcode, HugoSymbol symbol) {
        return table.contains(barcode, symbol);
    }

    /**
     * Returns the total number of records in this table.
     *
     * @return the total number of records in this table.
     */
    public int count() {
        return count;
    }

    /**
     * Counts the total number of missense-chop records for a given
     * tumor.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @return the total number of missense-chop records for the
     * specified tumor.
     */
    public int count(TumorBarcode barcode) {
        int total = 0;

        for (HugoSymbol symbol : table.viewInnerKeys(barcode))
            total += count(barcode, symbol);

        return total;
    }

    /**
     * Counts the number of missense-chop records for a given tumor
     * and gene.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @param symbol the gene of interest.
     *
     * @return the number of missense-chop records for the specified
     * tumor and gene.
     */
    public int count(TumorBarcode barcode, HugoSymbol symbol) {
        return lookup(barcode, symbol).size();
    }

    /**
     * Returns all missense-chop records for a given tumor.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @return an immutable list containing all missense-chop records
     * for the specified tumor (an empty list if there are no matching
     * records).
     */
    public List<MissCleavageRecord> lookup(TumorBarcode barcode) {
        List<MissCleavageRecord> records =
            new ArrayList<MissCleavageRecord>();

        for (HugoSymbol symbol : viewSymbols(barcode))
            records.addAll(lookup(barcode, symbol));

        return records;
    }

    /**
     * Returns all missense-chop records for a given tumor and gene.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @param symbol the HUGO symbol of interest.
     *
     * @return an immutable list containing all missense-chop records
     * for the specified tumor and gene (or an empty list if there are
     * no matching records).
     */
    public List<MissCleavageRecord> lookup(TumorBarcode barcode, HugoSymbol symbol) {
        RecordList recordList = table.get(barcode, symbol);

        if (recordList != null)
            return Collections.unmodifiableList(recordList);
        else
            return Collections.emptyList();
    }

    /**
     * Returns a read-only view of all tumor barcodes in this table.
     *
     * @return a read-only view of all tumor barcodes in this table.
     */
    public Set<TumorBarcode> viewBarcodes() {
        return table.viewOuterKeys();
    }

    /**
     * Returns a read-only view of all mutated genes for a given tumor.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @return a read-only view of all mutated genes for the specified
     * tumor.
     */
    public Set<HugoSymbol> viewSymbols(TumorBarcode barcode) {
        return table.viewInnerKeys(barcode);
    }
}
