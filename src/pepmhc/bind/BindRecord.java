
package pepmhc.bind;

import jean.peptide.Peptide;

/**
 * Provides the base class for all binding records managed by the
 * {@code pepmhc} project.
 */
public abstract class BindRecord {
    /**
     * The peptide described by this record.
     */
    protected final Peptide peptide;

    /**
     * Creates a new record describing a given peptide.
     *
     * @param peptide the peptided described by the new record.
     */
    protected BindRecord(Peptide peptide) {
        this.peptide = peptide;
    }

    /**
     * Returns the peptide described by this record.
     *
     * @return the peptide described by this record.
     */
    public final Peptide getPeptide() {
        return peptide;
    }
}
