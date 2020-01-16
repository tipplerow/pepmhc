
package pepmhc.calc;

import java.util.Collection;
import java.util.HashSet;

import jam.hla.Allele;
import jam.hla.Genotype;
import jam.math.DoubleUtil;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.binder.BindingThreshold;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;

/**
 * Computes presentation rates for MHC alleles and genotypes.
 */
public final class PresentationRateCalculator {
    private final PredictionMethod method;
    private final BindingThreshold threshold;
    private final Collection<Peptide> peptides;

    private static PresentationRateCalculator global = null;

    private PresentationRateCalculator(PredictionMethod method,
                                       BindingThreshold threshold,
                                       Collection<Peptide> peptides) {
        this.method = method;
        this.peptides = peptides;
        this.threshold = threshold;
    }

    /**
     * Returns a presentation rate calculator with a fixed prediction
     * method, binding threshold, and peptide collection.
     *
     * @param method the desired affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @param peptides the peptide collection to examine.
     *
     * @return a presentation rate calculator for the specified
     * prediction method, binding threshold, and peptide collection.
     */
    public static PresentationRateCalculator instance(PredictionMethod method,
                                                      BindingThreshold threshold,
                                                      Collection<Peptide> peptides) {
        return new PresentationRateCalculator(method, threshold, peptides);
    }

    /**
     * Computes the fraction of peptides that bind to an allele
     * according to a given prediction method and binding threshold.
     *
     * @param method the affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @param peptides the peptide collection to examine.
     *
     * @param allele the allele to test.
     *
     * @return the fraction of peptides that bind to the specified
     * allele according to the given prediction method and binding
     * threshold.
     */
    public static double compute(PredictionMethod method,
                                 BindingThreshold threshold,
                                 Collection<Peptide> peptides,
                                 Allele allele) {
        return instance(method, threshold, peptides).compute(allele);
    }

    /**
     * Computes the fraction of peptides that bind to a genotype
     * according to a given prediction method and binding threshold.
     *
     * @param method the affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @param peptides the peptide collection to examine.
     *
     * @param genotype the genotype to test.
     *
     * @return the fraction of peptides that bind to the specified
     * genotype according to the given prediction method and binding
     * threshold.
     */
    public static double compute(PredictionMethod method,
                                 BindingThreshold threshold,
                                 Collection<Peptide> peptides,
                                 Genotype genotype) {
        return instance(method, threshold, peptides).compute(genotype);
    }

    /**
     * Computes the fraction of peptides that bind to a given allele.
     *
     * @param allele the allele to test.
     *
     * @return the fraction of peptides that bind to the specified
     * allele.
     */
    public double compute(Allele allele) {
        Collection<BindingRecord> records = AffinityCache.get(method, allele, peptides);

        int bound = threshold.countBinders(records);
        int total = peptides.size();

        return DoubleUtil.ratio(bound, total);
    }

    /**
     * Computes the fraction of peptides that bind to a given
     * genotype.
     *
     * @param genotype the genotype to test.
     *
     * @return the fraction of peptides that bind to the specified
     * genotype.
     */
    public double compute(Genotype genotype) {
        //
        // A peptide may bind to multiple alleles, so to avoid
        // over-counting, collect the binder peptides from each
        // allele in a single set.
        //
        Collection<Peptide> binders = new HashSet<Peptide>();

        for (Allele allele : genotype.viewUniqueAlleles())
            binders.addAll(threshold.getBinders(AffinityCache.get(method, allele, peptides)));
        
        int bound = binders.size();
        int total = peptides.size();

        return DoubleUtil.ratio(bound, total);
    }

    /**
     * Computes the <em>ideal binding fraction</em> for a genotype:
     * the binding fraction the genotype would have if there was no
     * overlap in the binding repertoires of the alleles.
     *
     * @param genotype the genotype to test.
     *
     * @return the ideal binding fraction for the genotype.
     */
    public double computeIdeal(Genotype genotype) {
        double ideal = 0.0;

        for (Allele allele : genotype.viewUniqueAlleles())
            ideal += compute(allele);

        return ideal;
    }
}
