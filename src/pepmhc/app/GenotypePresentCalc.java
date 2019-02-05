
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

import pepmhc.calc.PresentationStatCalculator;

public final class GenotypePresentCalc extends JamApp {
    private final String genotypeInputFile;
    private final String alleleReportFile;
    private final String patientReportFile;

    private final PresentationStatCalculator calculator;

    // All unique alleles from all genotypes...
    private final Set<Allele> alleles;

    private String patientKeyName;
    private LineReader genotypeInputReader;
    private PrintWriter alleleReportWriter;
    private PrintWriter patientReportWriter;

    private GenotypePresentCalc(String... propFiles) {
        super(propFiles);

        this.genotypeInputFile = resolveGenotypeInputFile();
        this.alleleReportFile  = resolveAlleleReportFile();
        this.patientReportFile = resolvePatientReportFile();

        this.alleles = new TreeSet<Allele>();
        this.calculator = PresentationStatCalculator.global();
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
     * Name of the system property that specifies the full path name
     * of the genotype input file.
     *
     * The genotype file must contain two comma-separated columns. The
     * first column must contain a key identifying the patient to whom
     * the genotype belongs; the second column must list the alleles
     * contained in the genotype separated by white space.
     */
    public static final String GENOTYPE_INPUT_FILE_PROPERTY = "pepmhc.app.genotypeInputFile";

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

    private static String resolveGenotypeInputFile() {
        return JamProperties.getRequired(GENOTYPE_INPUT_FILE_PROPERTY);
    }

    private static String resolveAlleleReportFile() {   
        return JamProperties.getRequired(ALLELE_REPORT_FILE_PROPERTY);
    }

    private static String resolvePatientReportFile() {   
        return JamProperties.getRequired(PATIENT_REPORT_FILE_PROPERTY);
    }

    private void run() {
        processAlleles();
        processGenotypes();
    }

    private void processAlleles() {
        readAlleles();

        alleleReportWriter = IOUtil.openWriter(alleleReportFile);
        alleleReportWriter.println("allele,presentRate.mean,presentRate.sterr");

        try {
            for (Allele allele : alleles)
                reportAllele(allele);
        }
        finally {
            IOUtil.close(alleleReportWriter);
        }
    }
  
    private void readAlleles() {
        LineReader reader = LineReader.open(genotypeInputFile);

        try {
            // Skip header line...
            reader.next();

            for (String line : reader)
                alleles.addAll(parseGenotype(line));
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
        return Genotype.parse(columns[1], GENOTYPE_ALLELE_DELIM);
    }

    private void reportAllele(Allele allele) {
        StatSummary summary = calculator.compute(allele);

        alleleReportWriter.println(String.format("%s,%.8f,%.8f",
                                                 allele,
                                                 summary.getMean(),
                                                 summary.getError()));
    }

    private void processGenotypes() {
        genotypeInputReader = LineReader.open(genotypeInputFile);
        patientReportWriter = IOUtil.openWriter(patientReportFile);

        try {
            readGenotypeHeader();
            writeGenotypeHeader();

            for (String genotypeLine : genotypeInputReader)
                processGenotype(genotypeLine);
        }
        finally {
            IOUtil.close(genotypeInputReader);
            IOUtil.close(patientReportWriter);
        }
    }

    private void readGenotypeHeader() {
        patientKeyName = parsePatientKey(genotypeInputReader.next());
    }

    private void writeGenotypeHeader() {
        LineBuilder builder = LineBuilder.csv();

        builder.append(patientKeyName);
        builder.append("presentRate.mean");
        builder.append("presentRate.sterr");
        builder.append("idealRate.mean");
        builder.append("idealRate.sterr");

        patientReportWriter.println(builder);
    }

    private void processGenotype(String line) {
        String patientKey = parsePatientKey(line);
        Genotype genotype = parseGenotype(line);

        processGenotype(patientKey, genotype);
    }

    private void processGenotype(String patientKey, Genotype genotype) {
        JamLogger.info("Processing patient [%s]...", patientKey);

        StatSummary rateSummary  = calculator.compute(genotype);
        StatSummary idealSummary = calculator.computeIdeal(genotype);

        writeReportData(patientKey, rateSummary, idealSummary);
    }

    private void writeReportData(String patientKey, StatSummary rateSummary, StatSummary idealSummary) {
        LineBuilder builder = LineBuilder.csv();
        builder.append(patientKey);

        builder.append(rateSummary.getMean(), "%.8f");
        builder.append(rateSummary.getError(), "%.8f");
        builder.append(idealSummary.getMean(), "%.8f");
        builder.append(idealSummary.getError(), "%.8f");
    
        patientReportWriter.println(builder);
        patientReportWriter.flush();
    }

    public static void main(String[] propFiles) {
        GenotypePresentCalc calculator = new GenotypePresentCalc(propFiles);
        calculator.run();
    }
}
