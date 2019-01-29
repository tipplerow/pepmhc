
package pepmhc.engine;

import java.util.Collection;
import java.util.List;

import jam.lang.JamException;
import jam.peptide.Peptide;

import pepmhc.engine.net.NetMHCPredictor;
import pepmhc.engine.net.NetMHCPanPredictor;
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

        case SMM:
            return SMMPredictor.INSTANCE;

        case SMM_PMBEC:
            return SMMPmbecPredictor.INSTANCE;

        default:
            throw JamException.runtime("Unknown prediction method: [%s].", method);
        }
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
    public static BindingRecord predict(PredictionMethod method, String allele, Peptide peptide) {
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
    public static List<BindingRecord> predict(PredictionMethod method, String allele, Collection<Peptide> peptides) {
        return instance(method).predict(allele, peptides);
    }

    /**
     * Returns the prediction method for this predictor.
     *
     * @return the prediction method for this predictor.
     */
    public abstract PredictionMethod getMethod();

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
    public default BindingRecord predict(String allele, Peptide peptide) {
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
    public abstract List<BindingRecord> predict(String allele, Collection<Peptide> peptides);
}
