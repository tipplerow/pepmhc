
package pepmhc.chop;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import java.util.List;

import jam.app.JamLogger;
import jam.io.IOUtil;
import jam.lang.JamException;
import jam.math.Probability;

import jene.fasta.FastaPeptideRecord;
import jene.peptide.Peptide;

/**
 * Runs the {@code netchop} executable for a single peptide.
 */
public final class NetChopRunner {
    private final Peptide peptide;

    private File peptideFile;
    private List<Probability> cleavageScores;

    private NetChopRunner(Peptide peptide) {
        this.peptide = peptide;
    }

    /**
     * Runs the {@code netchop} executable for a single peptide to
     * compute the cleavage probabilities at each site in the peptide.
     *
     * @param peptide the peptide to score.
     *
     * @return a list containing the cleavage probabilities at each
     * site in the input peptide.
     */
    public static List<Probability> score(Peptide peptide) {
        NetChopRunner runner = new NetChopRunner(peptide);
        return runner.score();
    }

    private List<Probability> score() {
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
        FastaPeptideRecord record = new FastaPeptideRecord(peptideKey(), fastaComment(), peptide);

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
