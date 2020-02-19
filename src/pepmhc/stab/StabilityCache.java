
package pepmhc.stab;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jam.app.JamEnv;
import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.hla.Allele;
import jam.io.FileUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;
import jam.util.MapUtil;

/**
 * Maintains an in-memory and persistent cache of peptide-MHC
 * stability records indexed by HLA allele and peptide.
 */
public final class StabilityCache {
    private final Allele allele;
    private final Connection connection;
    private final Map<Peptide, StabilityRecord> records =
        new HashMap<Peptide, StabilityRecord>();

    private static final Map<Allele, StabilityCache> instances =
        new HashMap<Allele, StabilityCache>();

    private static final String TABLE_NAME = "stability";

    private static final String KEY_NAME = "peptide";
    private static final String KEY_TYPE = "string";

    private static final String HALF_LIFE_NAME = "half_life";
    private static final String HALF_LIFE_TYPE = "double";

    private static final String PERCENTILE_NAME = "percentile";
    private static final String PERCENTILE_TYPE = "double";

    private static final String TABLE_SCHEMA =
        String.format("%s %s PRIMARY KEY, %s %s, %s %s",
                      KEY_NAME, KEY_TYPE,
                      HALF_LIFE_NAME, HALF_LIFE_TYPE,
                      PERCENTILE_NAME, PERCENTILE_TYPE);

    private static final String CREATE_TABLE_COMMAND =
        String.format("CREATE TABLE IF NOT EXISTS %s (%s)", TABLE_NAME, TABLE_SCHEMA);

    private static final String SELECT_ALL_QUERY =
        String.format("SELECT %s, %s, %s FROM %s", KEY_NAME, HALF_LIFE_NAME, PERCENTILE_NAME, TABLE_NAME);

    private static final String INSERT_RECORD_TEMPLATE =
        "INSERT INTO %s VALUES('%s', %s, %s)";

    private StabilityCache(Allele allele) {
        this.allele = allele;
        this.connection = openConnection();

        loadTable();
    }

    private Connection openConnection() {
        try {
            Class.forName("org.sqlite.JDBC");

            Connection con = DriverManager.getConnection(formatURL());
            con.setAutoCommit(false);

            return con;
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
        String dirName;

        if (JamProperties.isSet(CACHE_DIRECTORY_PROPERTY))
            dirName = JamProperties.getRequired(CACHE_DIRECTORY_PROPERTY);
        else
            dirName = JamEnv.getRequired(CACHE_DIRECTORY_ENV);

        FileUtil.ensureDir(dirName);
        return dirName;
    }

    private String formatDbFile() {
        return formatAllele() + ".db";
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

        JamLogger.info("Loaded [%d] records from the stability database...", records.size());
    }

    private void loadRow(ResultSet resultSet) throws SQLException {
        Peptide peptide    = Peptide.parse(resultSet.getString(KEY_NAME));
        double  halfLife   = resultSet.getDouble(HALF_LIFE_NAME);
        double  percentile = resultSet.getDouble(PERCENTILE_NAME);

        cacheRecord(new StabilityRecord(peptide, halfLife, percentile));
    }

    private void cacheRecord(StabilityRecord record) {
        records.put(record.getPeptide(), record);
    }

    /**
     * Name of the environment variable that specifies the directory
     * containing the persistent database store. The system property
     * {@code pepmhc.stab.cacheDir} will take precedence if both are
     * specified.
     */
    public static final String CACHE_DIRECTORY_ENV = "PEPMHC_STABILITY_CACHE";

    /**
     * Name of the system property that specifies the directory
     * containing the persistent database store.
     */
    public static final String CACHE_DIRECTORY_PROPERTY = "pepmhc.stab.cacheDir";

    /**
     * Retrieves a binding record from the cache, computing on-demand
     * if necessary.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @param peptide the peptide bound to the MHC molecule.
     *
     * @return a stability record describing the peptide-MHC
     * interaction.
     */
    public static StabilityRecord get(Allele allele, Peptide peptide) {
        return get(allele, List.of(peptide)).get(0);
    }

    /**
     * Retrieves binding records from the cache, computing on-demand
     * if necessary.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @param peptides the peptides bound to the MHC molecule.
     *
     * @return the stability records describing the peptide-MHC
     * interactions.
     */
    public static List<StabilityRecord> get(Allele allele, Collection<Peptide> peptides) {
        return instance(allele).get(peptides);
    }

    private static StabilityCache instance(Allele allele) {
        JamLogger.info("Requesting stability cache for allele [%s]...", allele);
        StabilityCache cache = instanceSync(allele);

        JamLogger.info("Acquired stability cache for allele [%s]...", allele);
        return cache;
    }

    private static synchronized StabilityCache instanceSync(Allele allele) {
        StabilityCache cache = instances.get(allele);

        if (cache == null) {
            JamLogger.info("Creating stability cache for allele [%s]...", allele);
            cache = new StabilityCache(allele);
            instances.put(allele, cache);
        }

        return cache;
    }

    private List<StabilityRecord> get(Collection<Peptide> peptides) {
        JamLogger.info("Requesting [%d] records for allele [%s]...", peptides.size(), allele);
        List<StabilityRecord> records = getSync(peptides);

        JamLogger.info("Received [%d] records for allele [%s]...", peptides.size(), allele);
        return records;
    }

    private synchronized List<StabilityRecord> getSync(Collection<Peptide> peptides) {
        //
        // Identify peptides from the input collection that are not
        // present in the cache ("missing" peptides)...
        //
        Set<Peptide> missing = flagMissing(peptides);

        if (!missing.isEmpty()) {
            // Compute stability for each missing peptide...
            List<StabilityRecord> records = NetStab.run(allele, missing);

            // Add the new records to the cache and the database
            // table...
            updateCache(records);
            updateTable(records);
        }

        // Pull the records from the cache in the order requested...
        return MapUtil.get(records, peptides);
    }

    private Set<Peptide> flagMissing(Collection<Peptide> peptides) {
        Set<Peptide> missing = new HashSet<Peptide>();

        for (Peptide peptide : peptides)
            if (!records.containsKey(peptide))
                missing.add(peptide);

        return missing;
    }

    private void updateCache(List<StabilityRecord> records) {
        JamLogger.info("Adding [%s] records to the in-memory map...", records.size());

        for (StabilityRecord record : records)
            cacheRecord(record);
    }

    private void updateTable(List<StabilityRecord> records) {
        JamLogger.info("Adding [%s] records to the database table...", records.size());

        try (Statement statement = connection.createStatement()) {
            for (StabilityRecord record : records)
                insertRecord(statement, record);

            connection.commit();
        }
        catch (SQLException ex1) {
            //
            // No need to halt the program...
            //
            JamLogger.warn("Failed to update database table!");
            JamLogger.warn(ex1);

            try {
                connection.rollback();
            }
            catch (SQLException ex2) {
                JamLogger.warn("Failed to rollback update transaction!");
                JamLogger.warn(ex2);
            }
        }
    }

    private void insertRecord(Statement statement, StabilityRecord record) throws SQLException {
        String insertCommand =
            String.format(INSERT_RECORD_TEMPLATE,
                          TABLE_NAME,
                          formatPeptide(record.getPeptide()),
                          formatDouble(record.getHalfLife()),
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
