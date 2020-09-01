
package pepmhc.stab.proxy;

import java.util.Collection;
import java.util.List;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.affy.AffinityMethod;
import pepmhc.stab.StabilityMethod;
import pepmhc.stab.StabilityPredictor;
import pepmhc.stab.StabilityRecord;

/**
 * Predicts the stability of peptide-MHC complexes by computing
 * binding affinity with {@code netMHCpan} and using a regression
 * model to convert from affinity to half-life.
 */
public final class AffinityProxyPredictor extends StabilityPredictor {
    private AffinityProxyPredictor() {
    }

    /**
     * The single affinity-proxy predictor.
     */
    public static final AffinityProxyPredictor INSTANCE = new AffinityProxyPredictor();

    /**
     * The affinity method used to fit the proxy model.
     */
    public static final AffinityMethod AFFINITY_METHOD = AffinityMethod.NET_MHC_PAN;

    @Override public StabilityMethod getMethod() {
        return StabilityMethod.AFFINITY_PROXY;
    }

    @Override public boolean isInstalled() {
        //
        // The persistent parameter database is part of the
        // distribution, so there are models available...
        //
        return true;
    }

    @Override public List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return AffinityProxyStore.get(allele, AFFINITY_METHOD).predict(peptides);
    }
}
