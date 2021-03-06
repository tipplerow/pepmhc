
package pepmhc.engine.net;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;

import jam.app.JamLogger;
import jam.hla.Allele;
import jam.io.IOUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;

/**
 * Base class wrapper around a {@code netMHC} or {@code netMHCpan}
 * command-line process.
 */
public abstract class NetRunner {
    private final Allele allele;
    private final Collection<Peptide> peptides;

    private File peptideFile;
    private List<BindingRecord> records;

    /**
     * Creates a new runner for an allele and target peptides.
     *
     * @param allele the string code for the binding MHC allele.
     *
     * @param peptides the peptide targets.
     */
    protected NetRunner(Allele allele, Collection<Peptide> peptides) {
        this.allele = allele;
        this.peptides = peptides;
    }

    /**
     * Generates the exact and complete command-line request to pass
     * to the underlying {@code ProcessBuilder}.
     *
     * @param allele the string code for the binding MHC allele.
     *
     * @param peptideFile the input file containing the target peptides.
     *
     * @return the exact and complete command-line request to pass
     * to the underlying {@code ProcessBuilder}.
     */
    protected abstract List<String> formatCommand(Allele allele, File peptideFile);

    /**
     * Parses the output written by the command-line program.
     *
     * @param reader a reader for the command-line output stream.
     *
     * @return a list containing the binding records generated by the
     * command-line program.
     */
    protected abstract List<BindingRecord> parseOutput(BufferedReader reader);

    /**
     * Executes the command-line prediction process.
     *
     * @return the binding records generated by the command-line
     * prediction process.
     */
    protected List<BindingRecord> run() {
        JamLogger.info("Predicting binding affinity for [%d] peptides to allele [%s]...", peptides.size(), allele);

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

        return records;
    }

    private void initPeptideFile() throws IOException {
        peptideFile = File.createTempFile("NetRunner", ".pep");
    }
        
    private void writePeptideFile() throws IOException {
        PrintWriter writer = new PrintWriter(peptideFile);

        for (Peptide peptide : peptides)
            writer.println(peptide.formatString());

        writer.close();
    }

    private void launchProcess() throws IOException {
        ProcessBuilder builder = new ProcessBuilder(formatCommand(allele, peptideFile));
        Process        process = builder.start();
        BufferedReader reader  = IOUtil.openReader(process.getInputStream());

        records = parseOutput(reader);

        if (records.size() != peptides.size())
            throw JamException.runtime("Affinity prediction failed for allele [%s]!", allele);
    }
}
