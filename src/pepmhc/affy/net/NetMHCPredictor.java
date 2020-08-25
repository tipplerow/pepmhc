
package pepmhc.affy.net;

import java.io.File;
import java.util.Collection;
import java.util.List;

import jam.app.JamProperties;

import jean.hla.Allele;
import jean.peptide.Peptide;

import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityPredictor;
import pepmhc.affy.AffinityRecord;

public final class NetMHCPredictor extends AffinityPredictor {
    private NetMHCPredictor() {
    }

    /**
     * The single {@code netMHC} instance.
     */
    public static final NetMHCPredictor INSTANCE = new NetMHCPredictor();

    /**
     * Name of the environment variable that defines the full path to
     * the {@code netMHC} executable file.  If the system property
     * {@code pepmhc.netMHC} is also defined, it will take precedence.
     */
    public static final String EXECUTABLE_PATH_ENV = "NET_MHC_EXE";

    /**
     * Name of the system property that defines the full path to the
     * {@code netMHC} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.netMHC";

    /**
     * Resolves the full path to the {@code netMHC} executable file.
     *
     * @return the full path to the {@code netMHC} executable file.
     */
    public static File resolveExecutableFile() {
        return new File(resolveExecutableName());
    }

    /**
     * Resolves the full path to the {@code netMHC} executable file.
     *
     * @return the full path to the {@code netMHC} executable file.
     */
    public static String resolveExecutableName() {
        return JamProperties.resolve(EXECUTABLE_PATH_PROPERTY, EXECUTABLE_PATH_ENV, "netMHC");
    }

    @Override public AffinityMethod getMethod() {
        return AffinityMethod.NET_MHC;
    }

    @Override public boolean isInstalled() {
        return resolveExecutableFile().canExecute();
    }

    @Override public List<AffinityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return NetMHCRunner.run(allele, peptides);
    }
}
