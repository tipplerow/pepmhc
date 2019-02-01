
package pepmhc.engine.net;

import java.io.File;
import java.util.Collection;
import java.util.List;

import jam.app.JamProperties;
import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.engine.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class NetMHCPredictor implements Predictor {
    private NetMHCPredictor() {}

    /**
     * The single instance.
     */
    public static final NetMHCPredictor INSTANCE = new NetMHCPredictor();

    /**
     * Name of the system property that defines the full path to the
     * {@code netMHC} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.engine.net.netMHC";

    /**
     * Resolves the full path to the {@code netMHC} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY},
     * or {@code netMHC} in the working directory if the system property
     * is not set.
     */
    public static File resolveExecutableFile() {
        return new File(resolveExecutableName());
    }

    /**
     * Resolves the full path to the {@code netMHC} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY},
     * or {@code netMHC} in the working directory if the system property
     * is not set.
     */
    public static String resolveExecutableName() {
        return JamProperties.getRequired(EXECUTABLE_PATH_PROPERTY);
    }

    @Override public PredictionMethod getMethod() {
        return PredictionMethod.NET_MHC;
    }

    @Override public boolean isInstalled() {
        return resolveExecutableFile().canExecute();
    }

    @Override public List<BindingRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return NetMHCRunner.run(allele, peptides);
    }
}
