
package pepmhc.stab;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.lang.JamException;
import jam.peptide.Peptide;

import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class AffinityProxyPredictor extends StabilityPredictor {
    private final PredictionMethod affinityMethod;

    private AffinityProxyPredictor(PredictionMethod affinityMethod) {
        this.affinityMethod = affinityMethod;
    }

    /**
     * Returns the affinity-proxy prediction engine for a given
     * affinity prediction method.
     *
     * @param affinityMethod the affinity prediction method.
     *
     * @return the affinity-proxy prediction engine for the specified
     * affinity prediction method.
     */
    public static AffinityProxyPredictor instance(PredictionMethod affinityMethod) {
        return new AffinityProxyPredictor(affinityMethod);
    }

    @Override public StabilityMethod getMethod() {
        switch (affinityMethod) {
        case NET_MHC_PAN:
            return StabilityMethod.NET_MHC_PAN_AFFINITY_PROXY;

        default:
            throw JamException.runtime("Unknown affinity method: [%s].", affinityMethod);
        }
    }

    @Override public boolean isInstalled() {
        return Predictor.isInstalled(affinityMethod);
    }

    @Override public List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return getModel(allele).predict(peptides);
    }

    private AffinityProxyModel getModel(Allele allele) {
        return AffinityProxyModel.instance(allele, affinityMethod);
    }
}
