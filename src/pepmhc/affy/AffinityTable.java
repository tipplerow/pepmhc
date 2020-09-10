
package pepmhc.affy;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import jam.io.FileUtil;
import jam.math.Percentile;
import jam.sql.SQLColumn;
import jam.sql.SQLDb;
import jam.sql.SQLiteDb;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.bind.BindTable;

/**
 * Maintains a persistent database table of peptide-MHC affinity
 * records indexed by peptide.
 */
public final class AffinityTable extends BindTable<AffinityRecord> {
    private AffinityTable(SQLDb db) {
        super(db);
    }

    /**
     * The name of the {@code affinity} table.
     */
    public static final String TABLE_NAME = "affinity";

    /**
     * The name of the {@code affinity} column.
     */
    public static final String AFFINITY_NAME = "affinity";

    /**
     * Meta-data for the {@code affinity} column.
     */
    public static final SQLColumn AFFINITY_COLUMN =
        SQLColumn.create(AFFINITY_NAME, "double")
        .notNull();

    /**
     * Meta-data for the table columns.
     */
    public static final List<SQLColumn> COLUMN_LIST =
        List.of(PEPTIDE_COLUMN, AFFINITY_COLUMN, PERCENTILE_COLUMN);

    /**
     * Creates a new affinity table in an existing database.
     *
     * @param db the database that will contain the new table.
     *
     * @return a new affinity table contained in the specified
     * database.
     */
    public static AffinityTable create(SQLDb db) {
        return new AffinityTable(db);
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
    public static AffinityTable create(AffinityMethod method, Allele allele) {
        String dbFile = dbFile(method, allele);
        SQLDb  sqlDb  = SQLiteDb.instance(dbFile);

        return new AffinityTable(sqlDb);
    }

    /**
     * Returns the unique database file name for a given prediction
     * method and allele.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the unique database file name for the specified
     * prediction method and allele.
     */
    public static String dbFile(AffinityMethod method, Allele allele) {
        String cacheDir = AffinityCache.cacheDir();

        String dbDir = FileUtil.join(cacheDir, method.name());
        FileUtil.ensureDir(dbDir);

        return FileUtil.join(dbDir, formatAllele(allele) + ".db");
    }

    @Override public List<SQLColumn> getColumns() {
        return COLUMN_LIST;
    }

    @Override public AffinityRecord getRow(ResultSet resultSet) throws SQLException {
        Peptide peptide = Peptide.instance(resultSet.getString(1));
        Affinity affinity = Affinity.valueOf(getDouble(resultSet, 2));
        Percentile percentile = Percentile.valueOf(getDouble(resultSet, 3));

        return new AffinityRecord(peptide, affinity, percentile);
    }

    @Override public String getTableName() {
        return TABLE_NAME;
    }

    @Override public void prepareColumn(PreparedStatement statement, int index,
                                        AffinityRecord record, String columnName) throws SQLException {
        switch (columnName) {
        case PEPTIDE_NAME:
            statement.setString(index, record.getPeptide().formatString());
            break;

        case AFFINITY_NAME:
            statement.setDouble(index, record.getAffinity().doubleValue());
            break;

        case PERCENTILE_NAME:
            statement.setDouble(index, record.getPercentile().doubleValue());
            break;

        default:
            throw invalidColumn(columnName);
        }
    }
}
