
package pepmhc.cache;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import java.util.Collection;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.hla.Allele;
import jam.io.FileUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;
import jam.util.MapUtil;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

/**
 * Maintains an in-memory and persistent cache of peptide-MHC binding
 * affinities indexed by prediction method, HLA allele, and peptide.
 */
public final class AffinityCache {
    private final Allele allele;
    private final Predictor predictor;
    private final PredictionMethod method;
    private final Connection connection;
    private final Map<Peptide, BindingRecord> recordMap;

    private static final Map<PredictionMethod, Map<Allele, AffinityCache>> instances =
        new EnumMap<PredictionMethod, Map<Allele, AffinityCache>>(PredictionMethod.class);

    private static final String TABLE_NAME = "affinity";

    private static final String KEY_NAME = "peptide";
    private static final String KEY_TYPE = "string";

    private static final String AFFINITY_NAME = "affinity";
    private static final String AFFINITY_TYPE = "double";

    private static final String PERCENTILE_NAME = "percentile";
    private static final String PERCENTILE_TYPE = "double";

    private static final String TABLE_SCHEMA =
        String.format("%s %s PRIMARY KEY, %s %s, %s %s",
                      KEY_NAME, KEY_TYPE,
                      AFFINITY_NAME, AFFINITY_TYPE,
                      PERCENTILE_NAME, PERCENTILE_TYPE);

    private static final String CREATE_TABLE_COMMAND =
        String.format("CREATE TABLE IF NOT EXISTS %s (%s)", TABLE_NAME, TABLE_SCHEMA);

    private static final String SELECT_ALL_QUERY =
        String.format("SELECT %s, %s, %s FROM %s", KEY_NAME, AFFINITY_NAME, PERCENTILE_NAME, TABLE_NAME);

    private static final String INSERT_RECORD_TEMPLATE =
        "INSERT INTO %s VALUES('%s', %s, %s)";

    private AffinityCache(PredictionMethod method, Allele allele) {
        this.method = method;
        this.allele = allele;
        this.predictor = Predictor.instance(method);
        this.recordMap = new HashMap<Peptide, BindingRecord>();
        this.connection = openConnection();

        loadTable();
    }

    private Connection openConnection() {
        try {
            Class.forName("org.sqlite.JDBC");
            return DriverManager.getConnection(formatURL());
        }
        catch (Exception ex) {
            throw JamException.runtime(ex);
        }
    }

    private String formatURL() {
        return "jdbc:sqlite:" + resolveDbFile();
    }

    private String resolveDbFile() {
        return FileUtil.join(resolveDbDir(), formatDbFile());
    }

    private static String resolveDbDir() {
        return JamProperties.getRequired(CACHE_DIRECTORY_PROPERTY);
    }

    private String formatDbFile() {
        return String.format("%s_%s.db", method, formatAllele());
    }

    private String formatAllele() {
        return allele.longKey().replace("*", "-").replace(":", "-");
    }

    private void loadTable() {
        try (Statement statement = connection.createStatement()) {
            initTable(statement);
            loadTable(statement);
        }
        catch (SQLException ex) {
            throw JamException.runtime(ex);
        }
    }

    private void initTable(Statement statement) throws SQLException {
        //
        // Create the table if necessary...
        //
        statement.executeUpdate(CREATE_TABLE_COMMAND);
    }

    private void loadTable(Statement statement) throws SQLException {
        try (ResultSet resultSet = statement.executeQuery(SELECT_ALL_QUERY)) {
            while (resultSet.next()) {
                loadRow(resultSet);
            }
        }

        JamLogger.info("Loaded [%d] records from the affinity database...", recordMap.size());
    }

    private void loadRow(ResultSet resultSet) throws SQLException {
        Peptide peptide    = Peptide.parse(resultSet.getString(KEY_NAME));
        double  affinity   = resultSet.getDouble(AFFINITY_NAME);
        double  percentile = resultSet.getDouble(PERCENTILE_NAME);

        cacheRecord(new BindingRecord(peptide, affinity, percentile));
    }

    private void cacheRecord(BindingRecord record) {
        recordMap.put(record.getPeptide(), record);
    }

    /**
     * Name of the system property that specifies the directory
     * containing the persistent database store.
     */
    public static final String CACHE_DIRECTORY_PROPERTY = "pepmhc.cache.directory";

    /**
     * Retrieves a binding record from the cache, computing on-demand
     * if necessary.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @param peptide the peptide bound to the MHC molecule.
     *
     * @return a binding record describing the peptide-MHC
     * interaction.
     */
    public static BindingRecord get(PredictionMethod method, Allele allele, Peptide peptide) {
        return get(method, allele, List.of(peptide)).get(0);
    }

    /**
     * Retrieves binding records from the cache, computing on-demand
     * if necessary.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @param peptides the peptides bound to the MHC molecule.
     *
     * @return the binding records describing the peptide-MHC
     * interactions.
     */
    public static List<BindingRecord> get(PredictionMethod method, Allele allele, Collection<Peptide> peptides) {
        return instance(method, allele).get(peptides);
    }

    private static AffinityCache instance(PredictionMethod method, Allele allele) {
        Map<Allele, AffinityCache> cacheMap = instances.get(method);

        if (cacheMap == null) {
            JamLogger.info("Creating cache map for method [%s]...", method);
            cacheMap = new HashMap<Allele, AffinityCache>();
            instances.put(method, cacheMap);
        }

        AffinityCache cache = cacheMap.get(allele);

        if (cache == null) {
            JamLogger.info("Creating cache for allele [%s]...", allele);
            cache = new AffinityCache(method, allele);
            cacheMap.put(allele, cache);
        }

        return cache;
    }

    private List<BindingRecord> get(Collection<Peptide> peptides) {
        //
        // Identify peptides from the input collection that are not
        // present in the cache ("missing" peptides)...
        //
        Set<Peptide> missing = flagMissing(peptides);

        if (!missing.isEmpty()) {
            // Compute affinity for each missing peptide...
            List<BindingRecord> records = predictor.predict(allele, missing);

            // Add the new records to the cache and the database
            // table...
            updateCache(records);
            updateTable(records);
        }

        // Pull the records from the cache in the order requested...
        return MapUtil.get(recordMap, peptides);
    }

    private Set<Peptide> flagMissing(Collection<Peptide> peptides) {
        Set<Peptide> missing = new HashSet<Peptide>();

        for (Peptide peptide : peptides)
            if (!recordMap.containsKey(peptide))
                missing.add(peptide);

        return missing;
    }

    private void updateCache(List<BindingRecord> records) {
        JamLogger.info("Adding [%s] records to the in-memory map...", records.size());

        for (BindingRecord record : records)
            cacheRecord(record);
    }

    private void updateTable(List<BindingRecord> records) {
        JamLogger.info("Adding [%s] records to the database table...", records.size());

        try (Statement statement = connection.createStatement()) {
            for (BindingRecord record : records)
                insertRecord(statement, record);
        }
        catch (SQLException ex) {
            //
            // No need to halt the program...
            //
            JamLogger.warn("Failed to update database table!");
            JamLogger.warn(ex);
        }
    }

    private void insertRecord(Statement statement, BindingRecord record) throws SQLException {
        String insertCommand =
            String.format(INSERT_RECORD_TEMPLATE,
                          TABLE_NAME,
                          formatPeptide(record.getPeptide()),
                          formatDouble(record.getAffinity()),
                          formatDouble(record.getPercentile()));

        statement.executeUpdate(insertCommand);
    }

    private static String formatPeptide(Peptide peptide) {
        return peptide.formatString();
    }

    private static String formatDouble(double value) {
        if (Double.isNaN(value))
            return "null";
        else
            return Double.toString(value);
    }
}
