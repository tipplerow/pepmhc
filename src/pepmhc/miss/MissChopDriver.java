
package pepmhc.miss;

import java.util.List;

import jam.app.JamApp;
import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.io.IOUtil;

import jene.ensembl.EnsemblProteinDb;
import jene.hugo.HugoMaster;
import jene.missense.MissenseTable;
import jene.tcga.CellFraction;

/**
 * Processes a MAF file, generates the neo-peptides corresponding to
 * missense mutations, computes proteasomal cleavage probabilities for
 * the neo-peptides and their corresponding self-peptides, and writes
 * the records to an output file.
 */
public final class MissChopDriver extends JamApp {
    private final int peptideLength;
    private final String missChopFile;
    private final String hugoMasterFile;
    private final String missenseMAFFile;
    private final String ensemblProteomeFile;
    private final String ensemblSecondaryProteome;
    private final CellFraction ccfThreshold;

    private HugoMaster hugoMaster;
    private EnsemblProteinDb ensemblDb;

    private MissenseTable missenseTable;
    private List<MissChopRecord> missChopRecords;

    private MissChopDriver(String... propFiles) {
        super(propFiles);

        this.ccfThreshold = resolveCCFThreshold();
        this.missChopFile = resolveMissChopFile();
        this.peptideLength = resolvePeptideLength();
        this.hugoMasterFile = resolveHugoMasterFile();
        this.missenseMAFFile = resolveMissenseMAFFile();
        this.ensemblProteomeFile = resolveEnsemblProteomeFile();
        this.ensemblSecondaryProteome = resolveEnsemblSecondaryProteome();
    }

    private static int resolvePeptideLength() {
        return JamProperties.getRequiredInt(PEPTIDE_LENGTH_PROPERTY);
    }

    private static String resolveHugoMasterFile() {
        return JamProperties.getRequired(HUGO_MASTER_FILE);
    }

    private static String resolveMissChopFile() {
        return JamProperties.getRequired(MISS_CHOP_FILE_PROPERTY);
    }

    private static String resolveMissenseMAFFile() {
        return JamProperties.getRequired(MISSENSE_MAF_FILE_PROPERTY);
    }

    private static String resolveEnsemblProteomeFile() {
        return JamProperties.getRequired(ENSEMBL_PROTEOME_FILE);
    }

    private static String resolveEnsemblSecondaryProteome() {
        return JamProperties.getRequired(ENSEMBL_SECONDARY_PROTEOME);
    }

    private static CellFraction resolveCCFThreshold() {
        return CellFraction.valueOf(JamProperties.getRequired(CCF_THRESHOLD_PROPERTY));
    }

    /**
     * Name of the system property that specifies the minimum cancer
     * cell fraction required to process a missense mutation.
     */
    public static final String CCF_THRESHOLD_PROPERTY = "MissChopDriver.ccfThreshold";

    /**
     * Name of the system property that specifies the full path name
     * of the primary Ensembl protein database file.
     */
    public static final String ENSEMBL_PROTEOME_FILE = "MissChopDriver.ensemblProteomeFile";

    /**
     * Name of the system property that specifies the full path name
     * of the secondary Ensembl protein database file.
     */
    public static final String ENSEMBL_SECONDARY_PROTEOME = "MissChopDriver.ensemblSecondaryProteome";

    /**
     * Name of the system property that specifies the full path name
     * of the HUGO master file.
     */
    public static final String HUGO_MASTER_FILE = "MissChopDriver.hugoMasterFile";

    /**
     * Name of the system property that specifies the full path name
     * of the missense-chop output file.
     */
    public static final String MISS_CHOP_FILE_PROPERTY = "MissChopDriver.missChopFile";

    /**
     * Name of the system property that specifies the full path name
     * of the missense mutation MAF file.
     */
    public static final String MISSENSE_MAF_FILE_PROPERTY = "MissChopDriver.missenseMAFFile";

    /**
     * Name of the system property that specifies the length of the
     * peptide fragments to generate.
     */
    public static final String PEPTIDE_LENGTH_PROPERTY = "MissChopDriver.peptideLength";

    /**
     * Processes a MAF file, generates the neo-peptides corresponding to
     * missense mutations, computes proteasomal cleavage probabilities for
     * the neo-peptides and their corresponding self-peptides, and writes
     * the records to an output file.
     *
     * @param propFiles files containing the system properties that define
     * the runtime environment.
     *
     * @throws RuntimeException if any errors occur.
     */
    public static void run(String... propFiles) {
        MissChopDriver driver = new MissChopDriver(propFiles);
        driver.run();
    }

    private void run() {
        initializeEngine();
        loadMissenseTable();
        processMissenseTable();
        writeMissChopRecords();

        JamLogger.info("DONE!");
    }

    private void initializeEngine() {
        ensemblDb = EnsemblProteinDb.load(ensemblProteomeFile, ensemblSecondaryProteome);
        hugoMaster = HugoMaster.load(hugoMasterFile);

        MissChopEngine.initialize(hugoMaster, ensemblDb);
    }

    private void loadMissenseTable() {
        missenseTable = MissenseTable.load(missenseMAFFile, ccfThreshold);
    }

    private void processMissenseTable() {
        missChopRecords = MissChopEngine.generate(missenseTable, peptideLength);
    }

    private void writeMissChopRecords() {
        IOUtil.writeLines(missChopFile, false, MissChopRecord.header());
        IOUtil.writeObjects(missChopFile, true, missChopRecords, record -> record.format());
    }

    private static void usage() {
        System.err.println("Usage: jam.neo.MissChopDriver PROP_FILE1 [PROP_FILE2 ...]");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length < 1)
            usage();

        run(args);
    }
}
