
package pepmhc.stab;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.process.BatchProcess;

final class StabilityBatchProcess extends BatchProcess<Peptide, StabilityRecord> {
    private final Allele allele;

    private final static int BATCH_SIZE = 100000;

    private StabilityBatchProcess(Allele allele, Collection<Peptide> peptides) {
        super(peptides, BATCH_SIZE);
        this.allele = allele;
    }

    static List<StabilityRecord> run(Allele allele, Collection<Peptide> peptides) {
        StabilityBatchProcess process =
            new StabilityBatchProcess(allele, peptides);

        return process.runSequential();
    }

    @Override public List<StabilityRecord> runBatch(Collection<Peptide> inputSlice) {
        return NetStabRunner.run(allele, inputSlice);
    }
}

