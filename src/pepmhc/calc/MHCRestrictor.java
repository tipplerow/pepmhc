
package pepmhc.calc;

import java.util.Collection;
import java.util.HashSet;

import jam.hla.Allele;
import jam.hla.Genotype;
import jam.peptide.Peptide;
import jam.peptide.Peptidome;

import pepmhc.binder.BindingRecord;
import pepmhc.binder.BindingThreshold;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;

/**
 * Implements MHC restriction for HLA genotypes.
 */
public final class MHCRestrictor {
    private final PredictionMethod method;
    private final BindingThreshold threshold;

    private static MHCRestrictor global = null;

    private MHCRestrictor(PredictionMethod method, BindingThreshold threshold) {
        this.method = method;
        this.threshold = threshold;
    }

    /**
     * Returns the global MHC restrictor with the prediction method
     * and binding threshold specified by system properties.
     *
     * @return the global MHC restrictor.
     */
    public static MHCRestrictor global() {
        if (global == null)
            global = new MHCRestrictor(PredictionMethod.global(), BindingThreshold.global());

        return global;
    }

    /**
     * Returns an MHC restrictor with a fixed prediction method and
     * binding threshold.
     *
     * @param method the desired affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @return an MHC restrictor for the specified prediction method
     * and binding threshold.
     */
    public static MHCRestrictor instance(PredictionMethod method, BindingThreshold threshold) {
        return new MHCRestrictor(method, threshold);
    }

    /**
     * Implements the MHC restriction for a given genotype.
     *
     * @param genotype the HLA genotype of interest.
     *
     * @param peptides the set of peptides being restricted.
     *
     * @return all peptides that the genotype presents.
     */
    public Peptidome restrict(Genotype genotype, Collection<Peptide> peptides) {
        //
        // A peptide may bind to multiple alleles, so to avoid
        // over-counting, collect the binder peptides from each
        // allele in a single set.
        //
        Collection<Peptide> binders = new HashSet<Peptide>();

        for (Allele allele : genotype)
            binders.addAll(threshold.getBinders(AffinityCache.get(method, allele, peptides)));
        
        return Peptidome.create(binders);
    }
}
