
package pepmhc.stab;

import java.util.Collection;
import java.util.List;

import jam.app.JamProperties;
import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

/**
 * Enumerates peptide-MHC stability prediction methods.
 */
public enum StabilityMethod {
    /**
     * Direct prediction of half-life by the {@code netMHCstabpan}
     * engine.
     */
    NET_MHC_STAB_PAN {
        @Override public boolean isInstalled() {
            return NetStab.isInstalled();
        }

        @Override public List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides) {
            return StabilityCache.get(allele, peptides);
        }
    },

    /**
     * Inferred prediction of half-life by computing binding affinity
     * with {@code netMHCpan} and using a regression model to convert
     * from affinity to half-life.
     */
    NET_MHC_PAN_AFFINITY_PROXY {
        @Override public boolean isInstalled() {
            return Predictor.isInstalled(PredictionMethod.NET_MHC_PAN);
        }

        @Override public List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides) {
            return AffinityProxyModel.instance(allele, PredictionMethod.NET_MHC_PAN).predict(peptides);
        }
    };

    private static StabilityMethod global = null;

    /**
     * Name of the system property that specifies the global stability
     * prediction method.
     */
    public static final String METHOD_PROPERTY = "pepmhc.stab.method";

    /**
     * Returns the global prediction method specified through system
     * properties.
     *
     * @return the global prediction method.
     */
    public static StabilityMethod global() {
        if (global == null)
            global = JamProperties.getRequiredEnum(METHOD_PROPERTY, StabilityMethod.class);

        return global;
    }

    /**
     * Determines whether the underlying executable program or
     * prediction engine is available to the JVM.
     *
     * @return {@code true} iff the underlying executable program
     * or prediction engine is available to the JVM.
     */
    public abstract boolean isInstalled();

    /**
     * Predicts the stability for an allele and a single peptide.
     *
     * @param allele the code of the MHC allele presenting the
     * peptide.
     *
     * @param peptide the peptide being presented.
     *
     * @return the stability record for the allele and peptide.
     */
    public StabilityRecord predict(Allele allele, Peptide peptide) {
        return predict(allele, List.of(peptide)).get(0);
    }

    /**
     * Predicts the stability for an allele and a collection of
     * peptides.
     *
     * @param allele the code of the MHC allele presenting the
     * peptides.
     *
     * @param peptides the peptides being presented.
     *
     * @return a list of stability records for the allele and
     * peptides.
     */
    public abstract List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides);
}
