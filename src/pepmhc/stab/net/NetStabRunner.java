
package pepmhc.stab.net;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;

import jam.app.JamLogger;
import jam.io.IOUtil;
import jam.lang.JamException;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.stab.StabilityRecord;

/**
 * Executes {@code netMHCstabpan} command-line processes.
 */
public final class NetStabRunner {
    private final Allele allele;
    private final Collection<Peptide> peptides;

    private File peptideFile;
    private List<StabilityRecord> records;

    private NetStabRunner(Allele allele, Collection<Peptide> peptides) {
        this.allele = allele;
        this.peptides = peptides;
    }

    /**
     * Executes a {@code netMHCstabpan} command-line process.
     *
     * @param allele the binding MHC allele.
     *
     * @param peptides the peptide targets.
     *
     * @return the stability records generated by the command-line
     * prediction process.
     */
    public static List<StabilityRecord> run(Allele allele, Collection<Peptide> peptides) {
        NetStabRunner runner = new NetStabRunner(allele, peptides);
        return runner.run();
    }

    private List<StabilityRecord> run() {
        JamLogger.info("Predicting the stability of [%d] peptides for allele [%s]...", peptides.size(), allele);

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
        peptideFile = File.createTempFile("NetStabRunner", ".pep");
    }
        
    private void writePeptideFile() throws IOException {
        PrintWriter writer = new PrintWriter(peptideFile);

        for (Peptide peptide : peptides)
            writer.println(peptide.formatString());

        writer.close();
    }

    private void launchProcess() throws IOException {
        ProcessBuilder builder = new ProcessBuilder(formatCommand());
        Process        process = builder.start();
        BufferedReader reader  = IOUtil.openReader(process.getInputStream());

        try {
            records = NetStabParser.parse(reader);

            if (records.size() != peptides.size())
                throw JamException.runtime("Affinity prediction failed for allele [%s]!", allele);
        }
        finally {
            reader.close();
        }
    }

    private List<String>formatCommand() {
        return List.of(NetStab.resolveExecutableName(), 
                       "-a", formatAllele(),
                       "-p", peptideFile.getAbsolutePath());
    }

    private String formatAllele() {
        return allele.longKey().replace("*", "");
    }
}
