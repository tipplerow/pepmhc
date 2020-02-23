
package pepmhc.stab.proxy;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import jam.app.JamEnv;
import jam.hla.Allele;
import jam.io.FileUtil;
import jam.sql.SQLDb;
import jam.sql.SQLiteDb;
import jam.sql.SQLTable;

import pepmhc.affy.AffinityMethod;

/**
 * Maintains a parameter table for affinity-proxy models.
 */
public final class AffinityProxyTable extends SQLTable<AffinityProxyKey, AffinityProxyModel> {
    private AffinityProxyTable(SQLDb db) {
        super(db);
    }

    private static AffinityProxyTable instance;

    private static final String FILE_NAME  = "proxy_param.db";
    private static final String TABLE_NAME = "proxy_param";

    private static final String KEY_NAME = "allele_method";
    private static final String INTERCEPT_NAME = "intercept";
    private static final String COEFFICIENT_NAME = "coefficient";

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

    @Override public List<String> getColumnNames() {
        return List.of(KEY_NAME, INTERCEPT_NAME, COEFFICIENT_NAME);
    }

    @Override public AffinityProxyKey getKey(AffinityProxyModel model) {
        return model.getKey();
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

    @Override public String getTableSchema() {
        return String.format("%s string PRIMARY KEY, %s double, %s double",
                             KEY_NAME, INTERCEPT_NAME, COEFFICIENT_NAME);
    }

    @Override public void prepareInsert(PreparedStatement statement, AffinityProxyModel model) throws SQLException {
        statement.setString(1, model.getKey().format());
        statement.setDouble(2, model.getIntercept());
        statement.setDouble(3, model.getCoefficient());
    }
}
