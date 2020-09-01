
package pepmhc.stab.net;

import java.util.Collection;
import java.util.List;

import jam.process.BatchProcess;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.stab.StabilityRecord;

final class NetStabBatch extends BatchProcess<Peptide, StabilityRecord> {
    private final Allele allele;

    private final static int BATCH_SIZE = 100000;

    private NetStabBatch(Allele allele, Collection<Peptide> peptides) {
        super(peptides, BATCH_SIZE);
        this.allele = allele;
    }

    static List<StabilityRecord> run(Allele allele, Collection<Peptide> peptides) {
        NetStabBatch process =
            new NetStabBatch(allele, peptides);

        return process.runSequential();
    }

    @Override public List<StabilityRecord> runBatch(Collection<Peptide> inputSlice) {
        return NetStabRunner.run(allele, inputSlice);
    }
}

