
package pepmhc.miss;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamLogger;
import jam.io.LineReader;

import jene.tcga.TumorGeneRecordTable;

/**
 * Indexes missense-cleavage records by tumor barcode and HUGO symbol.
 */
public final class MissCleavageTable extends TumorGeneRecordTable<MissCleavageRecord> {
    private MissCleavageTable(Collection<MissCleavageRecord> records) {
        super(records);
    }

    /**
     * Populates a table from a collection of records.
     *
     * @param records the records to be indexed in the table.
     *
     * @return a table containing all records in the given collection.
     */
    public static MissCleavageTable create(Collection<MissCleavageRecord> records) {
        return new MissCleavageTable(records);
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
        return create(records);
    }
}
