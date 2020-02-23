
package pepmhc.engine.net;

import java.io.File;
import java.util.Collection;
import java.util.List;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class NetMHCPanPredictor implements Predictor {
    private NetMHCPanPredictor() {}

    /**
     * The single instance.
     */
    public static final NetMHCPanPredictor INSTANCE = new NetMHCPanPredictor();

    /**
     * Name of the environment variable that defines the full path to
     * the {@code netMHCpan} executable file.  If the system property
     * {@code pepmhc.engine.net.netMHCpan} is also defined, it will
     * take precedence.
     */
    public static final String EXECUTABLE_PATH_ENV = "NET_MHC_PAN_EXE";

    /**
     * Name of the system property that defines the full path to the
     * {@code netMHCpan} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.engine.net.netMHCpan";

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
        if (JamProperties.isSet(EXECUTABLE_PATH_PROPERTY))
            return JamProperties.getRequired(EXECUTABLE_PATH_PROPERTY);
        else
            return JamEnv.getRequired(EXECUTABLE_PATH_ENV);
    }

    @Override public PredictionMethod getMethod() {
        return PredictionMethod.NET_MHC_PAN;
    }

    @Override public boolean isInstalled() {
        return resolveExecutableFile().canExecute();
    }

    @Override public List<BindingRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return NetMHCPanRunner.run(allele, peptides);
    }
}
