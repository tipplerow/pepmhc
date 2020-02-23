
package pepmhc.cache;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.sql.SQLDb;
import jam.sql.SQLiteDb;
import jam.sql.SQLTable;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.PredictionMethod;

/**
 * Maintains a persistent database table of peptide-MHC affinity
 * records indexed by peptide.
 */
public final class AffinityTable extends SQLTable<Peptide, BindingRecord> {
    private AffinityTable(SQLDb db) {
        super(db);
    }

    private static final String TABLE_NAME = "affinity";

    private static final String KEY_NAME = "peptide";
    private static final String KEY_TYPE = "string";

    private static final String AFFINITY_NAME = "affinity";
    private static final String AFFINITY_TYPE = "double";

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
    public static AffinityTable create(PredictionMethod method, Allele allele) {
        String dbFile = AffinityDbFile.resolve(method, allele);
        SQLDb  sqlDb  = SQLiteDb.instance(dbFile);

        return new AffinityTable(sqlDb);
    }

    @Override public List<String> getColumnNames() {
        return List.of(KEY_NAME, AFFINITY_NAME, PERCENTILE_NAME);
    }

    @Override public Peptide getKey(BindingRecord record) {
        return record.getPeptide();
    }

    @Override public BindingRecord getRow(ResultSet resultSet) throws SQLException {
        Peptide peptide    = Peptide.parse(resultSet.getString(1));
        double  affinity   = getDouble(resultSet, 2);
        double  percentile = getDouble(resultSet, 3);

        return new BindingRecord(peptide, affinity, percentile);
    }

    @Override public String getTableName() {
        return TABLE_NAME;
    }

    @Override public String getTableSchema() {
        return String.format("%s %s PRIMARY KEY, %s %s, %s %s",
                             KEY_NAME, KEY_TYPE,
                             AFFINITY_NAME, AFFINITY_TYPE,
                             PERCENTILE_NAME, PERCENTILE_TYPE);
    }

    @Override public void prepareInsert(PreparedStatement statement, BindingRecord record) throws SQLException {
        statement.setString(1, record.getPeptide().formatString());
        statement.setDouble(2, record.getAffinity());
        statement.setDouble(3, record.getPercentile());
    }
}
