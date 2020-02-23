
package pepmhc.stab;

import pepmhc.bind.BindPredictor;

/**
 * Defines an interface for computational engines or algorithms that
 * predict the binding <em>stability</em> for peptide-MHC complexes:
 * the half-life for the dissociation of the complex.
 */
public abstract class StabilityPredictor extends BindPredictor<StabilityRecord> {
    /**
     * Returns the enumerated prediction method.
     *
     * @return the enumerated prediction method.
     */
    @Override public abstract StabilityMethod getMethod();
}
