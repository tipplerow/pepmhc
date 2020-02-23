
package pepmhc.affy;

import pepmhc.bind.BindPredictor;

/**
 * Defines an interface for computational engines or algorithms that
 * predict the binding <em>affinityy</em> for peptide-MHC complexes:
 * the IC50 concentration expressed in nanomolar units.
 */
public abstract class AffinityPredictor extends BindPredictor<AffinityRecord> {
    /**
     * Returns the enumerated prediction method.
     *
     * @return the enumerated prediction method.
     */
    @Override public abstract AffinityMethod getMethod();
}
