
package pepmhc.app;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.regex.Pattern;

import jam.app.JamApp;
import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.hla.Allele;
import jam.hla.Genotype;
import jam.io.IOUtil;
import jam.io.LineReader;
import jam.lang.JamException;
import jam.math.DoubleUtil;
import jam.math.StatSummary;
import jam.peptide.Peptide;
import jam.report.LineBuilder;
import jam.util.ListUtil;
import jam.util.RegexUtil;

import pepmhc.binder.BindingRecord;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;

public final class GenotypePresentCalc extends JamApp {
    private final int trialCount;
    private final int sampleSize;

    // Note that the PERCENTILE threshold takes precedence...
    private final double affinityThreshold;
    private final double percentileThreshold;

    private final String alleleReportFile;
    private final String patientReportFile;

    private final String peptideFileName;
    private final String genotypeFileName;

    private final PredictionMethod method;
    private final List<Set<Peptide>> peptideSamples;

    private final Set<Allele> alleleSet; // All individual alleles
    private final Map<Allele, StatSummary> alleleSummaries;

    // binderCache.get(trialIndex).get(allele) contains all peptides
    // from sample "trialIndex" that bind to "allele" with an affinity
    // below the threshold...
    private final List<Map<Allele, Set<Peptide>>> binderCache;

    private String patientKeyName;
    private LineReader genotypeReader;
    private PrintWriter alleleReportWriter;
    private PrintWriter patientReportWriter;

    private GenotypePresentCalc(String... propFiles) {
        super(propFiles);

        this.trialCount = resolveTrialCount();
        this.sampleSize = resolveSampleSize();

        this.affinityThreshold = resolveAffinityThreshold();
        this.percentileThreshold = resolvePercentileThreshold();

        this.alleleReportFile = resolveAlleleReportFile();
        this.patientReportFile = resolvePatientReportFile();

        this.peptideFileName = resolvePeptideFileName();
        this.genotypeFileName = resolveGenotypeFileName();

        this.method = resolvePredictionMethod();
        this.peptideSamples = new ArrayList<Set<Peptide>>();

        this.alleleSet = new TreeSet<Allele>();
        this.binderCache = new ArrayList<Map<Allele, Set<Peptide>>>();
        this.alleleSummaries = new HashMap<Allele, StatSummary>();
    }

    /**
     * Column delimiter for the genotype input file.
     */
    public static final Pattern GENOTYPE_COLUMN_DELIM = RegexUtil.COMMA;

    /**
     * Allele delimiter for the genotype input file.
     */
    public static final Pattern GENOTYPE_ALLELE_DELIM = RegexUtil.MULTI_WHITE_SPACE;

    /**
     * Name of the system property that specifies the affinity
     * threshold for peptide-MHC binding.
     */
    public static final String AFFINITY_THRESHOLD_PROPERTY = "pepmhc.app.affinityThreshold";

    /**
     * Name of the system property that specifies the percentile rank
     * threshold for peptide-MHC binding.  (The percentile threshold
     * takes precedence over the affinity threshold.)
     */
    public static final String PERCENTILE_THRESHOLD_PROPERTY = "pepmhc.app.percentileThreshold";

    /**
     * Name of the system property that specifies the full path name
     * of the genotype file.
     *
     * The genotype file must contain two comma-separated columns. The
     * first column must contain a key identifying the patient to whom
     * the genotype belongs; the second column must list the alleles
     * contained in the genotype separated by white space.
     */
    public static final String GENOTYPE_FILE_PROPERTY = "pepmhc.app.genotypeFile";

    /**
     * Name of the system property that specifies the full path name
     * to a flat file containing the peptides to examine.
     */
    public static final String PEPTIDE_FILE_PROPERTY = "pepmhc.app.peptideFile";

    /**
     * Name of the system property that specifies the peptide-MHC
     * affinity prediction method.
     */
    public static final String PREDICTION_METHOD_PROPERTY = "pepmhc.app.predictionMethod";

    /**
     * Name of the system property that specifies the base name of the
     * allele report file (written to the report directory.)
     */
    public static final String ALLELE_REPORT_FILE_PROPERTY = "pepmhc.app.alleleReportFile";

    /**
     * Name of the system property that specifies the base name of the
     * patient report file (written to the report directory.)
     */
    public static final String PATIENT_REPORT_FILE_PROPERTY = "pepmhc.app.patientReportFile";

    /**
     * Name of the system property that specifies the number of
     * peptides to include in each binding trial.
     */
    public static final String SAMPLE_SIZE_PROPERTY = "pepmhc.app.sampleSize";

    /**
     * Name of the system property that specifies the number of
     * independent binding trials to conduct.
     */
    public static final String TRIAL_COUNT_PROPERTY = "pepmhc.app.trialCount";

    private static double resolveAffinityThreshold() {
        return JamProperties.getRequiredDouble(AFFINITY_THRESHOLD_PROPERTY);
    }

    private static String resolveGenotypeFileName() {
        return JamProperties.getRequired(GENOTYPE_FILE_PROPERTY);
    }

    private static String resolvePeptideFileName() {
        return JamProperties.getRequired(PEPTIDE_FILE_PROPERTY);
    }

    private static PredictionMethod resolvePredictionMethod() {
        return JamProperties.getRequiredEnum(PREDICTION_METHOD_PROPERTY, PredictionMethod.class);
    }

    private static String resolveAlleleReportFile() {   
        return JamProperties.getRequired(ALLELE_REPORT_FILE_PROPERTY);
    }

    private static String resolvePatientReportFile() {   
        return JamProperties.getRequired(PATIENT_REPORT_FILE_PROPERTY);
    }

    private static double resolvePercentileThreshold() {
        return JamProperties.getRequiredDouble(PERCENTILE_THRESHOLD_PROPERTY);
    }

    private static int resolveSampleSize() {
        return JamProperties.getRequiredInt(SAMPLE_SIZE_PROPERTY);
    }

    private static int resolveTrialCount() {
        return JamProperties.getRequiredInt(TRIAL_COUNT_PROPERTY);
    }

    private void run() {
        samplePeptides();
        readAlleleSet();
        findBinders();
        processAlleles();
        processGenotypes();
    }
  
    private void samplePeptides() {
        List<Peptide> allPeptides = readPeptideFile();

        while (peptideSamples.size() < trialCount)
            peptideSamples.add(samplePeptides(allPeptides));
    }

    private List<Peptide> readPeptideFile() {
        JamLogger.info("Reading peptide file...");

        LineReader reader = LineReader.open(peptideFileName);
        List<Peptide> peptides = new ArrayList<Peptide>();

        try {
            for (String line : reader)
                peptides.add(Peptide.parse(line));
        }
        finally {
            IOUtil.close(reader);
        }

        return peptides;
    }

    private Set<Peptide> samplePeptides(List<Peptide> allPeptides) {
        JamLogger.info("Sampling [%d] of [%d] peptides...", sampleSize, allPeptides.size());

        Set<Peptide> sampledPeptides = new HashSet<Peptide>(sampleSize);

        while (sampledPeptides.size() < sampleSize)
            sampledPeptides.add(ListUtil.select(allPeptides));

        return sampledPeptides;
    }

    private void readAlleleSet() {
        LineReader reader = LineReader.open(genotypeFileName);

        try {
            // Skip header line...
            reader.next();

            for (String line : reader)
                alleleSet.addAll(parseGenotype(line));
        }
        finally {
            IOUtil.close(reader);
        }
    }

    private String parsePatientKey(String line) {
        String[] columns = RegexUtil.split(GENOTYPE_COLUMN_DELIM, line, 2);
        String   genoKey = columns[0];

        return genoKey;
    }

    private Genotype parseGenotype(String line) {
        String[] columns = RegexUtil.split(GENOTYPE_COLUMN_DELIM, line, 2);
        return Genotype.parse(line, GENOTYPE_ALLELE_DELIM);
    }

    private void findBinders() {
        for (int trialIndex = 0; trialIndex < trialCount; ++trialIndex)
            findBinders(trialIndex);
    }

    private void findBinders(int trialIndex) {
        JamLogger.info("Testing peptide sample [%d]...", trialIndex);

        if (binderCache.size() != trialIndex)
            throw new IllegalStateException("Inconsistent cache size.");

        binderCache.add(new HashMap<Allele, Set<Peptide>>());

        for (Allele allele : alleleSet)
            findBinders(trialIndex, allele);
    }

    private void findBinders(int trialIndex, Allele allele) {
        Set<Peptide> binderSet = new HashSet<Peptide>();

        List<BindingRecord> records =
            AffinityCache.get(method, allele, peptideSamples.get(trialIndex));

        if (records.size() != peptideSamples.get(trialIndex).size())
            throw JamException.runtime("Affinity computation failed for allele [%s]!", allele);

        for (BindingRecord record : records)
            if (isBound(record))
                binderSet.add(record.getPeptide());

        binderCache.get(trialIndex).put(allele, binderSet);
    }

    private boolean isBound(BindingRecord record) {
        return record.getPercentile() <= percentileThreshold || record.getAffinity() <= affinityThreshold;
    }

    private void processAlleles() {
        for (Allele allele : alleleSet)
            alleleSummaries.put(allele, processAllele(allele));

        writeAlleleReport();
    }

    private StatSummary processAllele(Allele allele) {
        return computeRateSummary(Set.of(allele));
    }

    private StatSummary computeRateSummary(Set<Allele> alleleSet) {
        return StatSummary.compute(computePresentationRates(alleleSet));
    }

    private double[] computePresentationRates(Set<Allele> alleleSet) {
        double[] rates = new double[trialCount];

        for (int trialIndex = 0; trialIndex < trialCount; ++trialIndex)
            rates[trialIndex] =
                computePresentationRate(trialIndex, alleleSet);

        return rates;
    }

    private double computePresentationRate(int trialIndex, Set<Allele> alleleSet) {
        int totalCount = peptideSamples.get(trialIndex).size();
        int binderCount = findBinderPeptides(trialIndex, alleleSet).size();

        return DoubleUtil.ratio(binderCount, totalCount);
    }

    private Set<Peptide> findBinderPeptides(int trialIndex, Set<Allele> alleleSet) {
        Set<Peptide> genotypeBinders = new HashSet<Peptide>();

        for (Allele allele : alleleSet)
            genotypeBinders.addAll(binderCache.get(trialIndex).get(allele));

        return genotypeBinders;
    }

    private void writeAlleleReport() {
        alleleReportWriter = IOUtil.openWriter(alleleReportFile);
        alleleReportWriter.println("allele,presentRate.mean,presentRate.sterr");

        for (Allele allele : alleleSet)
            alleleReportWriter.println(String.format("%s,%.8f,%.8f",
                                                     allele,
                                                     alleleSummaries.get(allele).getMean(),
                                                     alleleSummaries.get(allele).getError()));

        IOUtil.close(alleleReportWriter);
    }

    private void processGenotypes() {
        genotypeReader = LineReader.open(genotypeFileName);
        patientReportWriter = IOUtil.openWriter(patientReportFile);

        try {
            readGenotypeHeader();
            writeReportHeader();

            for (String genotypeLine : genotypeReader)
                processGenotype(genotypeLine);
        }
        finally {
            IOUtil.close(genotypeReader);
            IOUtil.close(patientReportWriter);
        }
    }

    private void readGenotypeHeader() {
        String   line    = genotypeReader.next();
        String[] columns = RegexUtil.split(GENOTYPE_COLUMN_DELIM, line, 2);

        patientKeyName = columns[0];
    }

    private void writeReportHeader() {
        LineBuilder builder = LineBuilder.csv();

        builder.append(patientKeyName);
        builder.append("presentRate.mean");
        builder.append("presentRate.sterr");
        builder.append("idealRate.mean");

        patientReportWriter.println(builder);
    }

    private void processGenotype(String line) {
        String patientKey = parsePatientKey(line);
        Genotype genotype = parseGenotype(line);

        processGenotype(patientKey, genotype);
    }

    private void processGenotype(String patientKey, Genotype genotype) {
        if (!isCovered(genotype)) {
            JamLogger.info("Incomplete allele coverage for patient [%s], skipping...", patientKey);
            return;
        }

        JamLogger.info("Processing patient [%s]...", patientKey);

        double idealRates = computeIdealRate(genotype);
        StatSummary rateSummary = computeRateSummary(genotype.viewUniqueAlleles());

        writeReportData(patientKey, idealRates, rateSummary);
    }

    private boolean isCovered(Genotype genotype) {
        for (Allele allele : genotype.viewUniqueAlleles())
            if (!binderCache.get(0).containsKey(allele))
                return false;

        return true;
    }

    private double computeIdealRate(Genotype genotype) {
        double idealRate = 0.0;

        for (Allele allele : genotype.viewUniqueAlleles())
            idealRate += alleleSummaries.get(allele).getMean();

        return idealRate;
    }

    private void writeReportData(String patientKey, double idealRate, StatSummary rateSummary) {
        LineBuilder builder = LineBuilder.csv();
        builder.append(patientKey);

        builder.append(rateSummary.getMean(), "%.8f");
        builder.append(rateSummary.getError(), "%.8f");
        builder.append(idealRate, "%.8f");
    
        patientReportWriter.println(builder);
        patientReportWriter.flush();
    }

    public static void main(String[] propFiles) {
        GenotypePresentCalc calculator = new GenotypePresentCalc(propFiles);
        calculator.run();
    }
}
