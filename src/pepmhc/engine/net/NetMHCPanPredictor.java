
package pepmhc.engine.net;

import java.io.File;
import java.util.Collection;
import java.util.List;

import jam.app.JamProperties;
import jam.peptide.Peptide;

import pepmhc.engine.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class NetMHCPanPredictor implements Predictor {
    private NetMHCPanPredictor() {}

    /**
     * The single instance.
     */
    public static final NetMHCPanPredictor INSTANCE = new NetMHCPanPredictor();

    /**
     * Name of the system property that defines the full path to the
     * {@code netMHCpan} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.engine.net.netMHCpan";

    /**
     * Determines whether the {@code netMHCpan} executable is installed
     * at the location specified by the {@code EXECUTABLE_PATH_PROPERTY}.
     *
     * @return {@code true} iff an executable file is installed at the
     * location specified by the {@code EXECUTABLE_PATH_PROPERTY}.
     */
    public static boolean isInstalled() {
        return resolveExecutableFile().canExecute();
    }

    /**
     * Resolves the full path to the {@code netMHCpan} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY},
     * or {@code netMHCpan} in the working directory if the system property
     * is not set.
     */
    public static File resolveExecutableFile() {
        return new File(resolveExecutableName());
    }

    /**
     * Resolves the full path to the {@code netMHCpan} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY},
     * or {@code netMHCpan} in the working directory if the system property
     * is not set.
     */
    public static String resolveExecutableName() {
        return JamProperties.getRequired(EXECUTABLE_PATH_PROPERTY);
    }

    @Override public PredictionMethod getMethod() {
        return PredictionMethod.NET_MHC_PAN;
    }

    @Override public List<BindingRecord> predict(String allele, Collection<Peptide> peptides) {
        return NetMHCPanRunner.run(allele, peptides);
    }
}
