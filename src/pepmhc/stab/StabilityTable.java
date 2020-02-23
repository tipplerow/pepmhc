
package pepmhc.stab;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import jam.hla.Allele;
import jam.io.FileUtil;
import jam.peptide.Peptide;
import jam.sql.SQLDb;
import jam.sql.SQLiteDb;
import jam.sql.SQLTable;

/**
 * Maintains a persistent database table of peptide-MHC stability
 * records indexed by peptide.
 */
public final class StabilityTable extends SQLTable<Peptide, StabilityRecord> {
    private StabilityTable(SQLDb db) {
        super(db);
    }

    private static final String TABLE_NAME = "stability";

    private static final String KEY_NAME = "peptide";
    private static final String KEY_TYPE = "string";

    private static final String HALF_LIFE_NAME = "half_life";
    private static final String HALF_LIFE_TYPE = "double";

    private static final String PERCENTILE_NAME = "percentile";
    private static final String PERCENTILE_TYPE = "double";

    /**
     * Creates a new affinity table in an existing database.
     *
     * @param db the database that will contain the new table.
     *
     * @return a new affinity table contained in the specified
     * database.
     */
    public static StabilityTable create(SQLDb db) {
        return new StabilityTable(db);
    }

    /**
     * Creates a new affinity table with a dedicated database for
     * a fixed allele and prediction method.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return a new affinity table for the specified allele and
     * prediction method.
     */
    public static StabilityTable create(StabilityMethod method, Allele allele) {
        String dbFile = dbFile(method, allele);
        SQLDb  sqlDb  = SQLiteDb.instance(dbFile);

        return new StabilityTable(sqlDb);
    }

    /**
     * Returns the unique database file name for a given prediction
     * method and allele.
     *
     * @param method the stability prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the unique database file name for the specified
     * prediction method and allele.
     */
    public static String dbFile(StabilityMethod method, Allele allele) {
        String cacheDir = StabilityCache.cacheDir();

        String dbDir = FileUtil.join(cacheDir, method.name());
        FileUtil.ensureDir(dbDir);

        return FileUtil.join(dbDir, formatAllele(allele) + ".db");
    }

    private static String formatAllele(Allele allele) {
        return allele.longKey().replace("*", "-").replace(":", "-");
    }

    @Override public List<String> getColumnNames() {
        return List.of(KEY_NAME, HALF_LIFE_NAME, PERCENTILE_NAME);
    }

    @Override public Peptide getKey(StabilityRecord record) {
        return record.getPeptide();
    }

    @Override public StabilityRecord getRow(ResultSet resultSet) throws SQLException {
        Peptide peptide    = Peptide.parse(resultSet.getString(1));
        double  halfLife   = getDouble(resultSet, 2);
        double  percentile = getDouble(resultSet, 3);

        return new StabilityRecord(peptide, halfLife, percentile);
    }

    @Override public String getTableName() {
        return TABLE_NAME;
    }

    @Override public String getTableSchema() {
        return String.format("%s %s PRIMARY KEY, %s %s, %s %s",
                             KEY_NAME, KEY_TYPE,
                             HALF_LIFE_NAME, HALF_LIFE_TYPE,
                             PERCENTILE_NAME, PERCENTILE_TYPE);
    }

    @Override public void prepareInsert(PreparedStatement statement, StabilityRecord record) throws SQLException {
        statement.setString(1, record.getPeptide().formatString());
        statement.setDouble(2, record.getHalfLife());
        statement.setDouble(3, record.getPercentile());
    }
}
