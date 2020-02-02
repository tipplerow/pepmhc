
package pepmhc.stab;

import java.sql.Connection;
import java.sql.ResultSet;

import jam.app.JamEnv;
import jam.app.JamLogger;
import jam.hla.Allele;
import jam.io.FileUtil;
import jam.lang.JamException;
import jam.util.SQLUtil;

import pepmhc.engine.PredictionMethod;

/**
 * Maintains a parameter table for affinity-proxy models.
 */
public final class AffinityProxyDb {
    private final Connection connection;
    private static AffinityProxyDb instance = null;

    private static final String KEY_DELIM = "|";

    private static final String FILE_NAME  = "proxy_param.db";
    private static final String TABLE_NAME = "proxy_param";

    private static final String KEY_NAME         = "allele_method";
    private static final String INTERCEPT_NAME   = "intercept";
    private static final String COEFFICIENT_NAME = "coefficient";

    private static final String TABLE_SCHEMA =
        String.format("%s string PRIMARY KEY, %s double, %s double",
                      KEY_NAME, INTERCEPT_NAME, COEFFICIENT_NAME);

    private AffinityProxyDb() {
        this.connection = SQLUtil.sqlite(formatURL());
        initTable();
    }

    private static String formatURL() {
        return FileUtil.join(JamEnv.getRequired("PEPMHC_HOME"), "data", "proxy", FILE_NAME);
    }

    private void initTable() {
        SQLUtil.createTable(connection, TABLE_NAME, TABLE_SCHEMA);
    }

    private static String key(AffinityProxyModel model) {
        return key(model.getAllele(), model.getMethod());
    }

    private static String key(Allele allele, PredictionMethod method) {
        return allele.shortKey() + KEY_DELIM + method.name();
    }

    /**
     * Returns the single database instance, created on demand.
     *
     * @return the single database instance.
     */
    public static AffinityProxyDb instance() {
        if (instance == null)
            instance = new AffinityProxyDb();

        return instance;
    }

    /**
     * Adds an affinity-proxy model to the parameter database (and
     * overwrites any previous parameters).
     *
     * @param model the model to add.
     */
    public void add(AffinityProxyModel model) {
        if (contains(model))
            remove(model);

        String updateCommand =
            String.format("INSERT INTO %s VALUES('%s', %.6f, %.6f)",
                          TABLE_NAME, key(model), model.getIntercept(), model.getCoefficient());

        JamLogger.info(updateCommand);
        SQLUtil.executeUpdate(connection, updateCommand);
    }

    /**
     * Determines whether the parameter table contains a row for an
     * affinity-proxy model.
     *
     * @param model the model of interest.
     *
     * @return {@code true} iff the parameter table contains a row for
     * the specified model.
     */
    public boolean contains(AffinityProxyModel model) {
        return contains(model.getAllele(), model.getMethod());
    }

    /**
     * Determines whether the parameter table contains a row for a
     * given allele and affinity prediction method.
     *
     * @param allele the allele of interest.
     *
     * @param method the affinity prediction method being used.
     *
     * @return {@code true} iff the parameter table contains a row for
     * the specified allele and affinity prediction method.
     */
    public boolean contains(Allele allele, PredictionMethod method) {
        int count = SQLUtil.count(connection, TABLE_NAME, KEY_NAME, key(allele, method));

        if (count < 0)
            throw JamException.runtime("Negative count.");

        if (count == 0)
            return false;

        if (count == 1)
            return true;

        throw JamException.runtime("Duplicate parameter records: [%s, %s].", allele, method);
    }

    /**
     * Returns the affinity-proxy model for a given allele and
     * affinity prediction method.
     *
     * @param allele the allele of interest.
     *
     * @param method the affinity prediction method being used.
     *
     * @return the affinity-proxy model with the parameters indexed by
     * the specified allele and prediction method (or {@code null} if
     * no such model has been stored in the database).
     */
    public AffinityProxyModel lookup(Allele allele, PredictionMethod method) {
        String queryString =
            String.format("SELECT %s, %s FROM %s WHERE %s = '%s'",
                          INTERCEPT_NAME, COEFFICIENT_NAME, TABLE_NAME, KEY_NAME, key(allele, method));

        JamLogger.info(queryString);
        ResultSet resultSet = SQLUtil.executeQuery(connection, queryString);

        double intercept;
        double coefficient;

        try {
            if (!resultSet.next())
                return null;

            intercept = resultSet.getDouble(1);
            coefficient = resultSet.getDouble(2);

            // Impossible...
            if (resultSet.next())
                throw JamException.runtime("Found duplicate parameter sets: [%s].", key(allele, method));
        }
        catch (Exception ex) {
            throw JamException.runtime(ex);
        }
        finally {
            SQLUtil.close(resultSet);
        }

        return new AffinityProxyModel(allele, method, intercept, coefficient);
    }

    /**
     * Removes the parameters of an affinity-proxy model from the table.
     *
     * @param model the model to remove.
     */
    public void remove(AffinityProxyModel model) {
        remove(model.getAllele(), model.getMethod());
    }

    /**
     * Removes the parameters of an affinity-proxy model for a given
     * allele and affinity prediction method.
     *
     * @param allele the allele to remove.
     *
     * @param method the affinity prediction method to remove.
     */
    public void remove(Allele allele, PredictionMethod method) {
        String updateCommand =
            String.format("DELETE FROM %s WHERE %s = '%s'",
                          TABLE_NAME, KEY_NAME, key(allele, method));

        JamLogger.info(updateCommand);
        SQLUtil.executeUpdate(connection, updateCommand);
    }
}
