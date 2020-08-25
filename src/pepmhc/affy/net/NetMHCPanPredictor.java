
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

public final class NetMHCPanPredictor extends AffinityPredictor {
    private NetMHCPanPredictor() {
    }

    /**
     * The single {@code netMHCpan} predictor.
     */
    public static final NetMHCPanPredictor INSTANCE = new NetMHCPanPredictor();

    /**
     * Name of the environment variable that defines the full path to
     * the {@code netMHCpan} executable file.  If the system property
     * {@code pepmhc.netMHCpan} is also defined, it takes precedence.
     */
    public static final String EXECUTABLE_PATH_ENV = "NET_MHC_PAN_EXE";

    /**
     * Name of the system property that defines the full path to the
     * {@code netMHCpan} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.netMHCpan";

    /**
     * Resolves the full path to the {@code netMHCpan} executable file.
     *
     * @return the full path to the {@code netMHCpan} executable file.
     */
    public static File resolveExecutableFile() {
        return new File(resolveExecutableName());
    }

    /**
     * Resolves the full path to the {@code netMHCpan} executable file.
     *
     * @return the full path to the {@code netMHCpan} executable file.
     */
    public static String resolveExecutableName() {
        return JamProperties.resolve(EXECUTABLE_PATH_PROPERTY, EXECUTABLE_PATH_ENV, "netMHCpan");
    }

    @Override public AffinityMethod getMethod() {
        return AffinityMethod.NET_MHC_PAN;
    }

    @Override public boolean isInstalled() {
        return resolveExecutableFile().canExecute();
    }

    @Override public List<AffinityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return NetMHCPanRunner.run(allele, peptides);
    }
}
