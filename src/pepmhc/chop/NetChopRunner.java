
package pepmhc.chop;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import java.util.List;

import jam.app.JamLogger;
import jam.fasta.FastaRecord;
import jam.io.IOUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;

/**
 * Runs the {@code netchop} executable for a single peptide.
 */
public final class NetChopRunner {
    private final Peptide peptide;

    private File peptideFile;
    private List<Double> cleavageScores;

    private NetChopRunner(Peptide peptide) {
        this.peptide = peptide;
    }

    /**
     * Runs the {@code netchop} executable for a single peptide.
     *
     * @param peptide the peptide to chop.
     *
     * @return a list containing the zero-based indexes of the
     * predicted cleavage sites.
     */
    public static List<Double> chop(Peptide peptide) {
        NetChopRunner runner = new NetChopRunner(peptide);
        return runner.chop();
    }

    private List<Double> chop() {
        try {
            initPeptideFile();
            writePeptideFile();
            launchProcess();
        }
        catch (IOException ioex) {
            JamLogger.error(ioex);
            throw JamException.runtime(ioex);
        }
        finally {
            peptideFile.delete();
        }

        return cleavageScores;
    }

    private void initPeptideFile() throws IOException {
        peptideFile = File.createTempFile("NetChopRunner", ".pep");
    }
        
    private void writePeptideFile() throws IOException {
        PrintWriter writer = new PrintWriter(peptideFile);
        FastaRecord record = new FastaRecord(peptideKey(), fastaComment(), peptide);

        writer.println(record.format());
        writer.close();
    }

    private static String peptideKey() {
        return "NetChop";
    }

    private static String fastaComment() {
        return "";
    }

    private void launchProcess() throws IOException {
        ProcessBuilder builder = new ProcessBuilder(formatCommand());
        Process        process = builder.start();
        BufferedReader reader  = IOUtil.openReader(process.getInputStream());

        cleavageScores = NetChopParser.parse(reader);
    }

    private List<String> formatCommand() {
        return List.of(NetChop.resolveExecutableName(), peptideFile.getAbsolutePath());
    }
}
