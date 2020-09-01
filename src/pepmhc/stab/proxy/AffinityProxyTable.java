
package pepmhc.stab.proxy;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import jam.app.JamEnv;
import jam.io.FileUtil;
import jam.sql.SQLColumn;
import jam.sql.SQLDb;
import jam.sql.SQLiteDb;
import jam.sql.SQLKeyTable;

import jene.hla.Allele;

import pepmhc.affy.AffinityMethod;

/**
 * Maintains a parameter table for affinity-proxy models.
 */
public final class AffinityProxyTable extends SQLKeyTable<AffinityProxyKey, AffinityProxyModel> {
    private AffinityProxyTable(SQLDb db) {
        super(db);
    }

    private static AffinityProxyTable instance;

    private static final String FILE_NAME  = "proxy_param.db";
    private static final String TABLE_NAME = "proxy_param";

    private static final String KEY_NAME = "allele_method";
    private static final String INTERCEPT_NAME = "intercept";
    private static final String COEFFICIENT_NAME = "coefficient";

    private static final SQLColumn KEY_COLUMN =
        SQLColumn.create(KEY_NAME, "string")
        .primaryKey();

    private static final SQLColumn INTERCEPT_COLUMN =
        SQLColumn.create(INTERCEPT_NAME, "double")
        .notNull();

    private static final SQLColumn COEFFICIENT_COLUMN =
        SQLColumn.create(COEFFICIENT_NAME, "double")
        .notNull();

    private static final List<SQLColumn> COLUMN_LIST =
        List.of(KEY_COLUMN, INTERCEPT_COLUMN, COEFFICIENT_COLUMN);

    /**
     * Returns the full path name of the persistent database file.
     *
     * @return the full path name of the persistent database file.
     */
    public static String dbFile() {
        return FileUtil.join(JamEnv.getRequired("PEPMHC_HOME"), "data", "proxy", FILE_NAME);
    }

    /**
     * Returns the single table instance, created on demand.
     *
     * @return the single table instance.
     */
    public static AffinityProxyTable instance() {
        if (instance == null)
            instance = new AffinityProxyTable(SQLiteDb.instance(dbFile()));

        return instance;
    }

    @Override public List<SQLColumn> getColumns() {
        return COLUMN_LIST;
    }

    @Override public AffinityProxyKey getKey(AffinityProxyModel model) {
        return model.getKey();
    }

    @Override public AffinityProxyKey getKey(ResultSet resultSet, String columnLabel) throws SQLException {
        return AffinityProxyKey.parse(resultSet.getString(columnLabel));
    }

    @Override public AffinityProxyModel getRow(ResultSet resultSet) throws SQLException {
        AffinityProxyKey key = AffinityProxyKey.parse(resultSet.getString(1));

        double intercept   = getDouble(resultSet, 2);
        double coefficient = getDouble(resultSet, 3);

        return new AffinityProxyModel(key, intercept, coefficient);
    }

    @Override public String getTableName() {
        return TABLE_NAME;
    }

    @Override public void prepareColumn(PreparedStatement statement, int index,
                                        AffinityProxyModel record, String columnName) throws SQLException {
        switch (columnName) {
        case KEY_NAME:
            statement.setString(index, record.getKey().format());
            break;

        case INTERCEPT_NAME:
            statement.setDouble(index, record.getIntercept());
            break;

        case COEFFICIENT_NAME:
            statement.setDouble(index, record.getCoefficient());
            break;

        default:
            throw invalidColumn(columnName);
        }
    }

    @Override public void prepareKey(PreparedStatement statement, int index, AffinityProxyKey key) throws SQLException {
        statement.setString(index, key.format());
    }
}
