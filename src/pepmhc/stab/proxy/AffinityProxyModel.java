
package pepmhc.stab.proxy;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jene.chem.HalfLife;
import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.affy.Affinity;
import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityRecord;
import pepmhc.affy.AffinityStore;

import pepmhc.stab.StabilityRecord;

public final class AffinityProxyModel {
    private final AffinityProxyKey key;
    private final double intercept;
    private final double coefficient;

    /**
     * Creates a new affinity-proxy model.
     *
     * @param key the unique key for the model.
     *
     * @param intercept the regression intercept for the model.
     *
     * @param coefficient the regression coefficient (slope) for the
     * model.
     *
     * @throws IllegalArgumentException unless the coefficient is
     * negative.
     */
    public AffinityProxyModel(AffinityProxyKey key, double intercept, double coefficient) {
        validateCoefficient(coefficient);

        this.key = key;
        this.intercept = intercept;
        this.coefficient = coefficient;
    }

    private static void validateCoefficient(double coefficient) {
        if (coefficient >= 0.0)
            throw new IllegalArgumentException("Coefficient must be negative.");
    }

    /**
     * Creates a new affinity-proxy model.
     *
     * @param allele the HLA allele which the model describes.
     *
     * @param method the affinity prediction method used to build the
     * model.
     *
     * @param intercept the regression intercept for the model.
     *
     * @param coefficient the regression coefficient (slope) for the
     * model.
     *
     * @throws IllegalArgumentException unless the coefficient is
     * negative.
     */
    public AffinityProxyModel(Allele allele, AffinityMethod method, double intercept, double coefficient) {
        this(AffinityProxyKey.instance(allele, method), intercept, coefficient);
    }

    /**
     * Returns the unique key for this model.
     *
     * @return the unique key for this model.
     */
    public AffinityProxyKey getKey() {
        return key;
    }

    /**
     * Returns the HLA allele which this model describes.
     *
     * @return the HLA allele which this model describes.
     */
    public Allele getAllele() {
        return key.getAllele();
    }

    /**
     * Returns the affinity prediction method used to build this
     * model.
     *
     * @return the affinity prediction method used to build this
     * model.
     */
    public AffinityMethod getMethod() {
        return key.getMethod();
    }

    /**
     * Returns the regression intercept for this model.
     *
     * @return the regression intercept for this model.
     */
    public double getIntercept() {
        return intercept;
    }

    /**
     * Returns the regression coefficient (slope) for this model.
     *
     * @return the regression coefficient (slope) for this model.
     */
    public double getCoefficient() {
        return coefficient;
    }

    /**
     * Computes the half-life that corresponds to a given affinity
     * according to this proxy model.
     *
     * @param affinity the peptide-MHC binding affinity.
     *
     * @return the half-life corresonding to the specified affinity
     * according to this proxy model.
     *
     * @throws IllegalArgumentException unless the affinity is positive.
     */
    public HalfLife halfLife(Affinity affinity) {
        //
        // Recall that the regression model is:
        //
        //     log(halfLife) = intercept + coefficient * log(affinity)
        //
        return HalfLife.valueOf(Math.exp(intercept + coefficient * Math.log(affinity.doubleValue())));
    }

    /**
     * Predicts the stability for a collection of peptides.
     *
     * @param peptides the peptides of interest.
     *
     * @return a list containing the stability records for the
     * specified peptides, in the same order as the collection
     * iterator.
     */
    public List<StabilityRecord> predict(Collection<Peptide> peptides) {
        List<AffinityRecord> affinityRecords =
            AffinityStore.instance(getMethod(), getAllele()).get(peptides);

        List<StabilityRecord> stabilityRecords =
            new ArrayList<StabilityRecord>(peptides.size());

        for (AffinityRecord affinityRecord : affinityRecords)
            stabilityRecords.add(stabilityRecord(affinityRecord));

        return stabilityRecords;
    }

    /**
     * Creates the stability record that corresponds to an affinity
     * record by converting the binding affinity to a dissociation
     * half-life.
     *
     * @param affinityRecord the affinity record to convert.
     *
     * @return the stability record that corresponds to the given
     * affinity record.
     */
    public StabilityRecord stabilityRecord(AffinityRecord affinityRecord) {
        return new StabilityRecord(affinityRecord.getPeptide(), halfLife(affinityRecord.getAffinity()));
    }
}
