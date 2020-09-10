
package pepmhc.miss;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamLogger;
import jam.io.LineReader;

import jene.tcga.TumorGeneRecordTable;

/**
 * Indexes missense-chop records by tumor barcode and HUGO symbol.
 */
public final class MissAffinityTable extends TumorGeneRecordTable<MissAffinityRecord> {
    private MissAffinityTable(Collection<MissAffinityRecord> records) {
        super(records);
    }

    /**
     * Populates a table from a collection of records.
     *
     * @param records the records to be indexed in the table.
     *
     * @return a table containing all records in the given collection.
     */
    public static MissAffinityTable create(Collection<MissAffinityRecord> records) {
        return new MissAffinityTable(records);
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
    public static MissAffinityTable load(String fileName) {
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
    public static MissAffinityTable load(File file) {
        List<MissAffinityRecord> records =
            new ArrayList<MissAffinityRecord>();

        try (LineReader reader = LineReader.open(file)) {
            // Skip header line...
            reader.next();

            for (String line : reader)
                records.add(MissAffinityRecord.parse(line));
        }

        JamLogger.info("MissAffinityTable: Loaded [%d] records.", records.size());
        return create(records);
    }
}
