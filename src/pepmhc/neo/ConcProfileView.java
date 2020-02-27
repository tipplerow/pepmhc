
package pepmhc.neo;

import java.util.Set;

import jam.hla.PeptideSource;
import jam.peptide.PeptideConcentrationProfile;

/**
 * Stores neo-peptide and self-peptide concentration profiles together
 * in a single object with read-only access.
 */
public final class ConcProfileView {
    private final PeptideConcentrationProfile neoProfile;
    private final PeptideConcentrationProfile selfProfile;

    private ConcProfileView(PeptideConcentrationProfile neoProfile,
                            PeptideConcentrationProfile selfProfile) {
        this.neoProfile = neoProfile;
        this.selfProfile = selfProfile;
    }

    /**
     * The single empty profile view.
     */
    public static final ConcProfileView EMPTY =
        create(PeptideConcentrationProfile.EMPTY,
               PeptideConcentrationProfile.EMPTY);

    /**
     * Creates a new concentration profile view from its constituent
     * profiles.
     *
     * @param neoProfile the concentration profile of neo-peptides.
     *
     * @param selfProfile the concentration profile of self-peptides.
     *
     * @return the new concentration profile view.
     */
    public static ConcProfileView create(PeptideConcentrationProfile neoProfile,
                                         PeptideConcentrationProfile selfProfile) {
        return new ConcProfileView(neoProfile, selfProfile);
    }

    /**
     * Returns the concentration profile for a given peptide source
     * type.
     *
     * @param source the peptide type: "neo" or "self".
     *
     * @return the concentration profile for peptides of the specified
     * type.
     */
    public PeptideConcentrationProfile getProfile(PeptideSource source) {
        switch (source) {
        case NEO:
            return neoProfile;

        case SELF:
            return selfProfile;

        default:
            throw new IllegalStateException("Unsupported peptide source.");
        }
    }

    /**
     * Returns the neo-peptide concentration profile.
     *
     * @return the neo-peptide concentration profile.
     */
    public PeptideConcentrationProfile getNeoProfile() {
        return neoProfile;
    }

    /**
     * Returns the self-peptide concentration profile.
     *
     * @return the self-peptide concentration profile.
     */
    public PeptideConcentrationProfile getSelfProfile() {
        return selfProfile;
    }
}
