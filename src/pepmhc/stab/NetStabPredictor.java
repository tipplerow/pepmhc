
package pepmhc.stab;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;

public final class NetStabPredictor extends StabilityPredictor {
    private NetStabPredictor() {
    }

    /**
     * Returns the single prediction engine.
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
