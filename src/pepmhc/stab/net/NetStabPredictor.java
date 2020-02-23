
package pepmhc.stab.net;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.stab.StabilityMethod;
import pepmhc.stab.StabilityPredictor;
import pepmhc.stab.StabilityRecord;

/**
 * Predicts the stability of peptide-MHC complexes using the 
 * {@code netMHCstabpan} program.
 */
public final class NetStabPredictor extends StabilityPredictor {
    private NetStabPredictor() {
    }

    /**
     * The single {@code netMHCstabpan} predictor.
     */
    public static final NetStabPredictor INSTANCE = new NetStabPredictor();

    @Override public StabilityMethod getMethod() {
        return StabilityMethod.NET_MHC_STAB_PAN;
    }

    @Override public boolean isInstalled() {
        return NetStab.isInstalled();
    }

    @Override public List<StabilityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        return NetStab.run(allele, peptides);
    }
}
