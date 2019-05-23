
package pepmhc.engine;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.lang.JamException;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.net.NetMHCPredictor;
import pepmhc.engine.net.NetMHCPanPredictor;
import pepmhc.engine.net.NetStabPredictor;
import pepmhc.engine.smm.SMMPredictor;
import pepmhc.engine.smm.SMMPmbecPredictor;

public interface Predictor {
    /**
     * Returns a predictor for a given prediction method.
     *
     * @param method the prediction method.
     *
     * @return the predictor with the specified prediction method.
     */
    public static Predictor instance(PredictionMethod method) {
        switch (method) {
        case NET_MHC:
            return NetMHCPredictor.INSTANCE;

        case NET_MHC_PAN:
            return NetMHCPanPredictor.INSTANCE;

        case NET_MHC_STAB_PAN:
            return NetStabPredictor.INSTANCE;

        case SMM:
            return SMMPredictor.INSTANCE;

        case SMM_PMBEC:
            return SMMPmbecPredictor.INSTANCE;

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
    public static boolean isInstalled(PredictionMethod method) {
        return instance(method).isInstalled();
    }

    /**
     * Predicts the binding between an allele and a single peptide.
     *
     * @param method the enumerated prediction method.
     *
     * @param allele the code of the MHC allele presenting the
     * peptide.
     *
     * @param peptide the peptide being presented.
     *
     * @return the binding record for the allele and peptide.
     */
    public static BindingRecord predict(PredictionMethod method, Allele allele, Peptide peptide) {
        return instance(method).predict(allele, peptide);
     }

    /**
     * Predicts the binding between an allele and a collection of
     * peptides.
     *
     * @param method the enumerated prediction method.
     *
     * @param allele the code of the MHC allele presenting the
     * peptides.
     *
     * @param peptides the peptides being presented.
     *
     * @return a list of binding records for the allele and peptides.
     */
    public static List<BindingRecord> predict(PredictionMethod method, Allele allele, Collection<Peptide> peptides) {
        return instance(method).predict(allele, peptides);
    }

    /**
     * Returns the prediction method for this predictor.
     *
     * @return the prediction method for this predictor.
     */
    public abstract PredictionMethod getMethod();

    /**
     * Determines whether the underlying executable program or
     * prediction engine is available to the JVM.
     *
     * @return {@code true} iff the underlying executable program
     * or prediction engine is available to the JVM.
     */
    public abstract boolean isInstalled();

    /**
     * Predicts the binding between an allele and a single peptide.
     *
     * @param allele the code of the MHC allele presenting the
     * peptide.
     *
     * @param peptide the peptide being presented.
     *
     * @return the binding record for the allele and peptide.
     */
    public default BindingRecord predict(Allele allele, Peptide peptide) {
        return predict(allele, List.of(peptide)).get(0);
    }

    /**
     * Predicts the binding between an allele and a collection of
     * peptides.
     *
     * @param allele the code of the MHC allele presenting the
     * peptides.
     *
     * @param peptides the peptides being presented.
     *
     * @return a list of binding records for the allele and peptides.
     */
    public abstract List<BindingRecord> predict(Allele allele, Collection<Peptide> peptides);
}
