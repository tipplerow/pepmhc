
package pepmhc.engine.net;

import java.io.BufferedReader;
import java.io.File;
import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;

/**
 * Executes {@code netMHCpan} command-line processes.
 */
public final class NetMHCPanRunner extends NetRunner {
    private NetMHCPanRunner(Allele allele, Collection<Peptide> peptides) {
        super(allele, peptides);
    }

    /**
     * Executes a {@code netMHCpan} command-line process.
     *
     * @param allele the string code for the binding MHC allele.
     *
     * @param peptides the peptide targets.
     *
     * @return the binding records generated by the command-line
     * prediction process.
     */
    public static List<BindingRecord> run(Allele allele, Collection<Peptide> peptides) {
        NetMHCPanRunner runner = new NetMHCPanRunner(allele, peptides);
        return runner.run();
    }

    @Override protected List<String>formatCommand(Allele allele, File peptideFile) {
        return List.of(NetMHCPanPredictor.resolveExecutableName(), 
                       "-a", formatAllele(allele),
                       "-BA", "-p", peptideFile.getAbsolutePath());
    }

    private static String formatAllele(Allele allele) {
        return allele.longKey().replace("*", "");
    }

    @Override protected List<BindingRecord> parseOutput(BufferedReader reader) {
        return NetMHCParser.parse(reader);
    }
}
