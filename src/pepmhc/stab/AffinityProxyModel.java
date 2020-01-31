
package pepmhc.stab;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.util.PairKeyTable;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class AffinityProxyModel {
    private final Allele allele;
    private final PredictionMethod method;

    private final double intercept;
    private final double coefficient;

    private static final PairKeyTable<Allele, PredictionMethod, AffinityProxyModel> instances = PairKeyTable.hash();

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
    public AffinityProxyModel(Allele allele, PredictionMethod method, double intercept, double coefficient) {
        validateCoefficient(coefficient);

        this.allele = allele;
        this.method = method;

        this.intercept = intercept;
        this.coefficient = coefficient;
    }

    private static void validateCoefficient(double coefficient) {
        if (coefficient >= 0.0)
            throw new IllegalArgumentException("Coefficient must be negative.");
    }

    /**
     * Returns the affinity-proxy prediction engine for a given HLA
     * allele and the default affinity prediction method.
     *
     * @param allele the HLA allele of interest.
     *
     * @return the affinity-proxy prediction engine for the specified
     * allele and the default affinity prediction method.
     */
    public static AffinityProxyModel instance(Allele allele) {
        return instance(allele, PredictionMethod.NET_MHC_PAN);
    }

    /**
     * Returns the affinity-proxy prediction engine for a given HLA
     * allele and affinity prediction method.
     *
     * @param allele the HLA allele of interest.
     *
     * @param method the affinity prediction method.
     *
     * @return the affinity-proxy prediction engine for the specified
     * allele and affinity prediction method.
     */
    public static AffinityProxyModel instance(Allele allele, PredictionMethod method) {
        AffinityProxyModel model = instances.get(allele, method);

        if (model == null) {
            model = create(allele, method);
            instances.put(allele, method, model);
        }

        return model;
    }

    private static AffinityProxyModel create(Allele allele, PredictionMethod method) {
        if (exists(allele, method))
            return load(allele, method);

        AffinityProxyModel model = AffinityProxyBuilder.build(allele, method);
        model.store();

        return model;
    }

    private static boolean exists(Allele allele, PredictionMethod method) {
        return AffinityProxyDb.instance().contains(allele, method);
    }

    private static AffinityProxyModel load(Allele allele, PredictionMethod method) {
        return AffinityProxyDb.instance().lookup(allele, method);
    }

    private void store() {
        AffinityProxyDb.instance().add(this);
    }

    /**
     * Returns the HLA allele which this model describes.
     *
     * @return the HLA allele which this model describes.
     */
    public Allele getAllele() {
        return allele;
    }

    /**
     * Returns the affinity prediction method used to build this
     * model.
     *
     * @return the affinity prediction method used to build this
     * model.
     */
    public PredictionMethod getMethod() {
        return method;
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
    public double halfLife(double affinity) {
        //
        // Recall that the regression model is:
        //
        //     log(halfLife) = intercept + coefficient * log(affinity)
        //
        if (affinity <= 0.0)
            throw new IllegalArgumentException("Negative affinity.");

        return Math.exp(intercept + coefficient * Math.log(affinity));
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
        List<BindingRecord> bindingRecords = Predictor.predict(method, allele, peptides);
        List<StabilityRecord> stabilityRecords = new ArrayList<StabilityRecord>(peptides.size());

        for (BindingRecord bindingRecord : bindingRecords)
            stabilityRecords.add(stabilityRecord(bindingRecord));

        return stabilityRecords;
    }

    /**
     * Creates the stability record that corresponds to a binding
     * record by converting the binding affinity to a dissociation
     * half-life.
     *
     * @param bindingRecord the binding record to convert.
     *
     * @return the stability record that corresponds to the given
     * binding record.
     */
    public StabilityRecord stabilityRecord(BindingRecord bindingRecord) {
        return new StabilityRecord(bindingRecord.getPeptide(), halfLife(bindingRecord.getAffinity()));
    }
}
