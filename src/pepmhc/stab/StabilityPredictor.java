
package pepmhc.stab;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.lang.JamException;
import jam.peptide.Peptide;

import pepmhc.engine.PredictionMethod;

public abstract class StabilityPredictor {
    private static StabilityPredictor global = null;

    /**
     * Returns the predictor for the global prediction method.
     *
     * @return the predictor for the global prediction method.
     */
    public static StabilityPredictor global() {
        if (global == null)
            global = instance(StabilityMethod.global());

        return global;
    }

    /**
     * Returns a predictor for a given prediction method.
     *
     * @param method the prediction method.
     *
     * @return the predictor with the specified prediction method.
     */
    public static StabilityPredictor instance(StabilityMethod method) {
        switch (method) {
        case NET_MHC_STAB_PAN:
            return NetStabPredictor.INSTANCE;

        case NET_MHC_PAN_AFFINITY_PROXY:
            return AffinityProxyPredictor.instance(PredictionMethod.NET_MHC_PAN);

        default:
            throw JamException.runtime("Unknown prediction method: [%s].", method);
        }
    }

    /**
     * Determines whether a given prediction method is installed and
     * available for the JVM.
     *
     * @param method the prediction method.
     *
     * @return {@code true} iff the specified prediction method is
     * installed and available for the JVM.
     */
    public static boolean isInstalled(StabilityMethod method) {
        return instance(method).isInstalled();
    }

    /**
     * Predicts the stability for an allele and a single peptide.
     *
     * @param method the enumerated prediction method.
     *
     * @param allele the code of the MHC allele presenting the
     * peptide.
     *
     * @param peptide the peptide being presented.
     *
     * @return the stability record for the allele and peptide.
     */
    public static StabilityRecord predict(StabilityMethod method, Allele allele, Peptide peptide) {
        return instance(method).predict(allele, peptide);
     }

    /**
     * Predicts the stability for an allele and a collection of
     * peptides.
     *
     * @param method the enumerated prediction method.
     *
     * @param allele the code of the MHC allele presenting the
     * peptides.
     *
     * @param peptides the peptides being presented.
     *
     * @return a list of stability records for the allele and peptides.
     */
    public static List<StabilityRecord> predict(StabilityMethod method, Allele allele, Collection<Peptide> peptides) {
        return instance(method).predict(allele, peptides);
    }

    /**
     * Returns the prediction method for this predictor.
     *
     * @return the prediction method for this predictor.
     */
    public abstract StabilityMethod getMethod();

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
