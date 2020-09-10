
package pepmhc.miss;

import java.util.List;

import jam.app.JamApp;
import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.io.IOUtil;

import jene.tcga.TumorGenotypeTable;

import pepmhc.affy.AffinityMethod;

/**
 * Computes the MHC binding affinity for neo-peptide and self-peptide
 * pairs contained in cleavage records.
 */
public final class MissAffinityDriver extends JamApp {
    private final String missAffinityFile;
    private final String missCleavageFile;
    private final String tumorPatientFile;
    private final String patientGenotypeFile;
    private final AffinityMethod affinityMethod;

    private MissCleavageTable cleavageTable;
    private TumorGenotypeTable genotypeTable;
    private List<MissAffinityRecord> affinityRecords;

    private MissAffinityDriver(String... propFiles) {
        super(propFiles);

        this.affinityMethod = resolveAffinityMethod();
        this.missAffinityFile = resolveMissAffinityFile();
        this.missCleavageFile = resolveMissCleavageFile();
        this.tumorPatientFile = resolveTumorPatientFile();
        this.patientGenotypeFile = resolvePatientGenotypeFile();
    }

    private static AffinityMethod resolveAffinityMethod() {
        return JamProperties.getRequiredEnum(AFFINITY_METHOD_PROPERTY, AffinityMethod.class);
    }

    private static String resolveMissAffinityFile() {
        return JamProperties.getRequired(MISS_AFFINITY_FILE_PROPERTY);
    }

    private static String resolveMissCleavageFile() {
        return JamProperties.getRequired(MISS_CLEAVAGE_FILE_PROPERTY);
    }

    private static String resolveTumorPatientFile() {
        return JamProperties.getRequired(TUMOR_PATIENT_FILE_PROPERTY);
    }

    private static String resolvePatientGenotypeFile() {
        return JamProperties.getRequired(PATIENT_GENOTYPE_FILE_PROPERTY);
    }

    /**
     * Name of the system property that specifies the enumerated
     * affinity prediction method.
     */
    public static final String AFFINITY_METHOD_PROPERTY = "MissAffinityDriver.affinityMethod";

    /**
     * Name of the system property that specifies the full path name
     * of the affinity record output file.
     */
    public static final String MISS_AFFINITY_FILE_PROPERTY = "MissAffinityDriver.missAffinityFile";

    /**
     * Name of the system property that specifies the full path name
     * of the cleavage record input file.
     */
    public static final String MISS_CLEAVAGE_FILE_PROPERTY = "MissAffinityDriver.missCleavageFile";

    /**
     * Name of the system property that specifies the full path name
     * of the patient genotype file.
     */
    public static final String PATIENT_GENOTYPE_FILE_PROPERTY = "MissAffinityDriver.patientGenotypeFile";

    /**
     * Name of the system property that specifies the length of the
     * tumor patient mapping file.
     */
    public static final String TUMOR_PATIENT_FILE_PROPERTY = "MissAffinityDriver.tumorPatientFile";

    /**
     * Computes the MHC binding affinity for neo-peptide and self-peptide
     * pairs contained in cleavage records.
     *
     * @param propFiles files containing the system properties that define
     * the runtime environment.
     *
     * @throws RuntimeException if any errors occur.
     */
    public static void run(String... propFiles) {
        MissAffinityDriver driver = new MissAffinityDriver(propFiles);
        driver.run();
    }

    private void run() {
        loadCleavageTable();
        loadGenotypeTable();
        genAffinityRecords();
        writeAffinityRecords();

        JamLogger.info("DONE!");
    }

    private void loadCleavageTable() {
        cleavageTable = MissCleavageTable.load(missCleavageFile);
    }

    private void loadGenotypeTable() {
        genotypeTable = TumorGenotypeTable.load(tumorPatientFile, patientGenotypeFile);
    }

    private void genAffinityRecords() {
        affinityRecords =
            MissAffinityEngine.generate(affinityMethod,
                                        cleavageTable,
                                        genotypeTable);
    }

    private void writeAffinityRecords() {
        IOUtil.writeLines(missAffinityFile, false, MissAffinityRecord.header());
        IOUtil.writeObjects(missAffinityFile, true, affinityRecords, record -> record.format());
    }

    private static void usage() {
        System.err.println("Usage: jam.neo.MissAffinityDriver PROP_FILE1 [PROP_FILE2 ...]");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length < 1)
            usage();

        run(args);
    }
}
