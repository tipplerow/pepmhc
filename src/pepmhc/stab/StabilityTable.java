
package pepmhc.stab;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import jam.io.FileUtil;
import jam.sql.SQLColumn;
import jam.sql.SQLDb;
import jam.sql.SQLiteDb;

import jean.hla.Allele;
import jean.peptide.Peptide;

import pepmhc.bind.BindTable;

/**
 * Maintains a persistent database table of peptide-MHC stability
 * records indexed by peptide.
 */
public final class StabilityTable extends BindTable<StabilityRecord> {
    private StabilityTable(SQLDb db) {
        super(db);
    }

    /**
     * The name of the {@code stability} table.
     */
    public static final String TABLE_NAME = "stability";

    /**
     * The name of the {@code half_life} column.
     */
    public static final String HALF_LIFE_NAME = "half_life";

    /**
     * Meta-data for the {@code half_life} column.
     */
    public static final SQLColumn HALF_LIFE_COLUMN =
        SQLColumn.create(HALF_LIFE_NAME, "double")
        .notNull();

    /**
     * Meta-data for the table columns.
     */
    public static final List<SQLColumn> COLUMN_LIST =
        List.of(PEPTIDE_COLUMN, HALF_LIFE_COLUMN, PERCENTILE_COLUMN);

    /**
     * Creates a new stability table in an existing database.
     *
     * @param db the database that will contain the new table.
     *
     * @return a new stability table contained in the specified
     * database.
     */
    public static StabilityTable create(SQLDb db) {
        return new StabilityTable(db);
    }

    /**
     * Creates a new stability table with a dedicated database for
     * a fixed allele and prediction method.
     *
     * @param method the stability prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return a new stability table for the specified allele and
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

    @Override public List<SQLColumn> getColumns() {
        return COLUMN_LIST;
    }

    @Override public StabilityRecord getRow(ResultSet resultSet) throws SQLException {
        Peptide peptide    = Peptide.instance(resultSet.getString(1));
        double  halfLife   = getDouble(resultSet, 2);
        double  percentile = getDouble(resultSet, 3);

        return new StabilityRecord(peptide, halfLife, percentile);
    }

    @Override public String getTableName() {
        return TABLE_NAME;
    }

    @Override public void prepareColumn(PreparedStatement statement, int index,
                                        StabilityRecord record, String columnName) throws SQLException {
        switch (columnName) {
        case PEPTIDE_NAME:
            statement.setString(index, record.getPeptide().formatString());
            break;

        case HALF_LIFE_NAME:
            statement.setDouble(index, record.getHalfLife());
            break;

        case PERCENTILE_NAME:
            statement.setDouble(index, record.getPercentile());
            break;

        default:
            throw invalidColumn(columnName);
        }
    }
}
