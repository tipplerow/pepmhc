
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
import jam.util.SetUtil;

/**
 * Maintains a persistent database table of peptide-MHC stability
 * records indexed by HLA allele and peptide.
 */
public final class StabilityStore {
    private final Allele allele;

    private static final Map<Allele, StabilityStore> instances =
        new HashMap<Allele, StabilityStore>();

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

    private StabilityStore(Allele allele) {
        this.allele = allele;
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
     * @return a stability record describing the peptide-MHC
     * interaction.
     */
    public static synchronized StabilityStore instance(Allele allele) {
        StabilityStore instance = instances.get(allele);

        if (instance == null) {
            instance = new StabilityStore(allele);
            instances.put(allele, instance);
        }

        return instance;
    }

    /**
     * Returns the allele served by this stability store.
     *
     * @return the allele served by this stability store.
     */
    public Allele getAllele() {
        return allele;
    }

    /**
     * Retrieves a binding record from the store, computing it
     * on-demand if necessary.
     *
     * @param peptide the peptide bound to the MHC molecule.
     *
     * @return a stability record describing the peptide-MHC
     * interaction.
     */
    public StabilityRecord getRecord(Peptide peptide) {
        return getRecords(List.of(peptide)).get(0);
    }

    /**
     * Retrieves binding records from the store, computing them
     * on-demand if necessary.
     *
     * @param peptides the peptides bound to the MHC molecule.
     *
     * @return the stability records describing the peptide-MHC
     * interactions.
     */
    public List<StabilityRecord> getRecords(Collection<Peptide> peptides) {
        List<StabilityRecord> records;

        try {
            JamLogger.info("Requesting [%d] stability records for allele [%s]...", peptides.size(), allele.longKey());
            records = getRecordsSync(peptides);
        }
        catch (Exception ex) {
            throw JamException.runtime(ex);
        }

        JamLogger.info("Received [%d] stability records for allele [%s].", records.size(), allele.longKey());
        return records;
    }

    private synchronized List<StabilityRecord> getRecordsSync(Collection<Peptide> requestedPeptides)
        throws ClassNotFoundException, SQLException {
        //
        // Load all previously computed stability records...
        //
        Map<Peptide, StabilityRecord> tableRecords = loadTable();

        // Identify peptides from the input collection that are not
        // present in the table...
        Set<Peptide> missingPeptides = flagMissing(tableRecords.keySet(), requestedPeptides);

        if (!missingPeptides.isEmpty()) {
            // Compute stability for each missing peptide...
            List<StabilityRecord> computedRecords = NetStab.run(allele, missingPeptides);

            // Add the new records to the database table and our local map...
            updateMap(tableRecords, computedRecords);
            updateTable(computedRecords);
        }

        return MapUtil.get(tableRecords, requestedPeptides);
    }

    private Map<Peptide, StabilityRecord> loadTable() throws ClassNotFoundException, SQLException {
        try (Connection connection = openConnection()) {
            return loadTable(connection);
        }
    }

    private Connection openConnection() throws ClassNotFoundException, SQLException {
        Class.forName("org.sqlite.JDBC");

        Connection con = DriverManager.getConnection(formatURL());
        con.setAutoCommit(false);

        return con;
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

    private Map<Peptide, StabilityRecord> loadTable(Connection connection) throws SQLException {
        Map<Peptide, StabilityRecord> records;
        JamLogger.info("Loading the [%s] stability database...", allele.longKey());

        try (Statement statement = connection.createStatement()) {
            initTable(statement);
            records = loadTable(statement);
        }

        JamLogger.info("Loaded [%d] records from the [%s] stability database.", records.size(), allele.longKey());
        return records;
    }

    private static void initTable(Statement statement) throws SQLException {
        //
        // Create the table if necessary...
        //
        statement.executeUpdate(CREATE_TABLE_COMMAND);
    }

    private static Map<Peptide, StabilityRecord> loadTable(Statement statement) throws SQLException {
        Map<Peptide, StabilityRecord> records =
            new HashMap<Peptide, StabilityRecord>();

        try (ResultSet resultSet = statement.executeQuery(SELECT_ALL_QUERY)) {
            while (resultSet.next()) {
                StabilityRecord record = loadRow(resultSet);
                records.put(record.getPeptide(), record);
            }
        }

        return records;
    }

    private static StabilityRecord loadRow(ResultSet resultSet) throws SQLException {
        Peptide peptide    = Peptide.parse(resultSet.getString(KEY_NAME));
        double  halfLife   = resultSet.getDouble(HALF_LIFE_NAME);
        double  percentile = resultSet.getDouble(PERCENTILE_NAME);

        return new StabilityRecord(peptide, halfLife, percentile);
    }

    private static Set<Peptide> flagMissing(Set<Peptide> tablePeptides, Collection<Peptide> requestedPeptides) {
        Set<Peptide> missingPeptides = new HashSet<Peptide>();

        for (Peptide peptide : requestedPeptides)
            if (!tablePeptides.contains(peptide))
                missingPeptides.add(peptide);

        return missingPeptides;
    }

    private static void updateMap(Map<Peptide, StabilityRecord> tableRecords, List<StabilityRecord> computedRecords) {
        for (StabilityRecord computedRecord : computedRecords)
            tableRecords.put(computedRecord.getPeptide(), computedRecord);
    }

    private void updateTable(List<StabilityRecord> records) throws ClassNotFoundException, SQLException {
        JamLogger.info("Adding [%s] records to the [%s] stability database...", records.size(), allele.longKey());

        try (Connection connection = openConnection()) {
            updateTable(connection, records);
        }

        JamLogger.info("Added [%s] records to the [%s] stability database.", records.size(), allele.longKey());
    }

    private static void updateTable(Connection connection, List<StabilityRecord> records) {
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

    private static void insertRecord(Statement statement, StabilityRecord record) throws SQLException {
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
