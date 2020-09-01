
package pepmhc.calc;

import java.util.ArrayList;
import java.util.List;

import jam.math.StatSummary;

import jene.hla.Allele;
import jene.hla.Genotype;

import pepmhc.affy.AffinityThreshold;
import pepmhc.affy.AffinityMethod;

/**
 * Computes presentation rate statistics for MHC alleles and
 * genotypes taken over a collection of peptide samples.
 */
public final class PresentationStatCalculator {
    private final List<PresentationRateCalculator> rateCalcs;

    private static PresentationStatCalculator global = null;

    private PresentationStatCalculator(AffinityMethod method,
                                       AffinityThreshold threshold,
                                       PeptideSample    peptides) {
        int sampleCount = peptides.countSamples();

        this.rateCalcs =
            new ArrayList<PresentationRateCalculator>(sampleCount);

        for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex)
            rateCalcs.add(PresentationRateCalculator.instance(method, threshold, peptides.viewSample(sampleIndex)));
    }

    /**
     * Returns a presentation statistic calculator for the global
     * prediction method, binding threshold, and peptide sample.
     *
     * @return a presentation statistic calculator for the global
     * prediction method, binding threshold, and peptide sample.
     */
    public static PresentationStatCalculator global() {
        if (global == null)
            global = new PresentationStatCalculator(AffinityMethod.global(),
                                                    AffinityThreshold.global(),
                                                    PeptideSample.global());
        return global;
    }

    /**
     * Returns a presentation statistic calculator with a fixed
     * prediction method, binding threshold, and peptide sample set.
     *
     * @param method the desired affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @param peptides a reference peptide set to examine.
     *
     * @return a presentation statistic calculator for the specified
     * prediction method, binding threshold, and peptide set.
     */
    public static PresentationStatCalculator instance(AffinityMethod method,
                                                      AffinityThreshold threshold,
                                                      PeptideSample    peptides) {
        return new PresentationStatCalculator(method, threshold, peptides);
    }

    /**
     * Computes a statistical summary for the fraction of peptides
     * that bind to an allele according to a prediction method and
     * binding threshold.
     *
     * @param method the affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @param peptides the reference peptide collections to examine.
     *
     * @param allele the allele to test.
     *
     * @return a summary for the the fraction of peptides that bind to
     * the specified allele according to the given prediction method
     * and binding threshold.
     */
    public static StatSummary compute(AffinityMethod method,
                                      AffinityThreshold threshold,
                                      PeptideSample peptides,
                                      Allele allele) {
        return instance(method, threshold, peptides).compute(allele);
    }

    /**
     * Computes a statistical summary for the fraction of peptides
     * that bind to a genotype according to a prediction method and
     * binding threshold.
     *
     * @param method the affinity prediction method.
     *
     * @param threshold the affinity/percentile threshold for binding.
     *
     * @param peptides the reference peptide collections to examine.
     *
     * @param genotype the genotype to test.
     *
     * @return a summary for the fraction of peptides that bind to the
     * specified genotype according to the given prediction method and
     * binding threshold.
     */
    public static StatSummary compute(AffinityMethod method,
                                      AffinityThreshold threshold,
                                      PeptideSample peptides,
                                      Genotype genotype) {
        return instance(method, threshold, peptides).compute(genotype);
    }

    /**
     * Computes a statistical summary for the fraction of reference
     * peptides that bind to a given allele.
     *
     * @param allele the allele to test.
     *
     * @return a statistical summary for the fraction of reference
     * peptides that bind to the specified allele.
     */
    public StatSummary compute(Allele allele) {
        List<Double> rates = new ArrayList<Double>();

        for (PresentationRateCalculator rateCalc : rateCalcs)
            rates.add(rateCalc.compute(allele));

        return StatSummary.compute(rates);
    }

    /**
     * Computes a statistical summary for the fraction of reference
     * peptides that bind to a given genotype.
     *
     * @param genotype the genotype to test.
     *
     * @return a statistical summary for the fraction of reference
     * peptides that bind to the specified genotype.
     */
    public StatSummary compute(Genotype genotype) {
        List<Double> rates = new ArrayList<Double>();

        for (PresentationRateCalculator rateCalc : rateCalcs)
            rates.add(rateCalc.compute(genotype));

        return StatSummary.compute(rates);
    }

    /**
     * Computes a statistical summary for the ideal binding fraction
     * of a genotype: the binding fraction the genotype would have if
     * there was no overlap in the binding repertoires of the alleles.
     *
     * @param genotype the genotype to test.
     *
     * @return a statistical summary of the ideal binding fraction for
     * the specified genotype.
     */
    public StatSummary computeIdeal(Genotype genotype) {
        List<Double> rates = new ArrayList<Double>();

        for (PresentationRateCalculator rateCalc : rateCalcs)
            rates.add(rateCalc.computeIdeal(genotype));

        return StatSummary.compute(rates);
    }
}
