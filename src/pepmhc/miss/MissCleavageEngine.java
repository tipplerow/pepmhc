
package pepmhc.miss;

import java.util.ArrayList;
import java.util.List;

import jam.app.JamLogger;
import jam.math.UnitIndex;
import jam.math.UnitIndexRange;
import jam.math.Probability;
import jam.util.ListUtil;
import jam.util.StreamUtil;

import jene.ensembl.EnsemblProteinDb;
import jene.hugo.HugoMaster;
import jene.hugo.HugoSymbol;
import jene.missense.MissenseGroup;
import jene.missense.MissenseRecord;
import jene.missense.MissenseTable;
import jene.neo.NeoPeptide;
import jene.neo.SelfPeptide;
import jene.peptide.Peptide;
import jene.peptide.ProteinChange;
import jene.tcga.TumorBarcode;

import pepmhc.chop.NetChopEngine;

/**
 * Generates the self/neo-peptide pairs corresponding to missense
 * mutations and computes their proteasomal cleavage probabilities.
 */
public final class MissCleavageEngine {
    private final int peptideLength;
    private final HugoSymbol hugoSymbol;
    private final TumorBarcode tumorBarcode;
    private final MissenseGroup missenseGroup;

    private Peptide nativeProtein;
    private Peptide mutatedProtein;

    private NetChopEngine nativeChopEngine;
    private NetChopEngine mutatedChopEngine;

    private static HugoMaster hugoMaster = null;
    private static EnsemblProteinDb ensemblDb = null;

    private MissCleavageEngine(MissenseGroup missenseGroup, int peptideLength) {
        this.missenseGroup = missenseGroup;
        this.peptideLength = peptideLength;

        this.hugoSymbol = missenseGroup.getHugoSymbol();
        this.tumorBarcode = missenseGroup.getTumorBarcode();
    }

    /**
     * Assigns the shared data structures that are used to process all
     * mutation groups.  This method must be called before any groups
     * are processed.
     *
     * @param hugoMaster the mapping from HUGO symbols to Ensembl genes.
     *
     * @param ensemblDb the Ensembl protein database.
     */
    public static void initialize(HugoMaster hugoMaster, EnsemblProteinDb ensemblDb) {
        MissCleavageEngine.ensemblDb = ensemblDb;
        MissCleavageEngine.hugoMaster = hugoMaster;
    }

    /**
     * Generates the self/neo-peptide pairs corresponding to a group
     * of missense mutations and computes their proteasomal cleavage
     * probabilities.
     *
     * @param missenseGroup a group of missense mutations observed
     * in the same tumor sample and gene.
     *
     * @param peptideLength the desired length of the self-peptide
     * and neo-peptide fragments.
     *
     * @return a list of cleavage records for the input mutation
     * group.
     *
     * @throws RuntimeException if the Ensembl database and HUGO
     * master have not been initialized or if the native peptide
     * cannot be resolved.
     */
    public static List<MissCleavageRecord> generate(MissenseGroup missenseGroup, int peptideLength) {
        if (!isInitialized())
            throw new IllegalStateException("The MissCleavageEngine has not been initialized.");

        MissCleavageEngine engine =
            new MissCleavageEngine(missenseGroup, peptideLength);

        try {
            return engine.generate();
        }
        catch (RuntimeException ex) {
            JamLogger.warn(ex);
            return List.of();
        }
    }

    /**
     * Generates the self/neo-peptide pairs corresponding to each
     * group of missense mutations in a patient cohort and computes
     * their proteasomal cleavage probabilities.
     *
     * @param missenseTable a table of missense mutations observed
     * in a patient cohort.
     *
     * @param peptideLength the desired length of the self-peptide
     * and neo-peptide fragments.
     *
     * @return a list of cleavage records for the input mutation
     * table.
     *
     * @throws RuntimeException if the Ensembl database and HUGO
     * master have not been initialized or if any native peptides
     * cannot be resolved.
     */
    public static List<MissCleavageRecord> generate(MissenseTable missenseTable, int peptideLength) {
        List<MissenseGroup> missenseGroups = missenseTable.group();

        List<List<MissCleavageRecord>> groupRecords =
            StreamUtil.applyParallel(missenseGroups, group -> generate(group, peptideLength));

        List<MissCleavageRecord> missCleavageRecords =
            ListUtil.cat(groupRecords);

        JamLogger.info("Sorting cleavage records...");
        missCleavageRecords.sort(MissCleavageRecord.COMPARATOR);

        return missCleavageRecords;
    }

    private static boolean isInitialized() {
        return ensemblDb != null && hugoMaster != null;
    }

    private List<MissCleavageRecord> generate() {
        JamLogger.info("Generating cleavage records: [%s, %s]...",
                       tumorBarcode.getKey(), hugoSymbol.getKey());

        nativeProtein = missenseGroup.resolveNative(ensemblDb, hugoMaster);
        mutatedProtein = missenseGroup.mutate(nativeProtein);

        assert mutatedProtein.length() == nativeProtein.length();

        nativeChopEngine = NetChopEngine.run(nativeProtein);
        mutatedChopEngine = NetChopEngine.run(mutatedProtein);

        List<MissCleavageRecord> missCleavageRecords =
            new ArrayList<MissCleavageRecord>();

        for (MissenseRecord missenseRecord : missenseGroup)
            missCleavageRecords.addAll(generate(missenseRecord));

        return missCleavageRecords;
    }

    private List<MissCleavageRecord> generate(MissenseRecord missenseRecord) {
        ProteinChange proteinChange = missenseRecord.getProteinChange();
        int proteinChangePosition = proteinChange.getPosition().getUnitIndex();

        List<MissCleavageRecord> missCleavageRecords =
            new ArrayList<MissCleavageRecord>(peptideLength);
        
        for (int neoPepMissPos = 1; neoPepMissPos <= peptideLength; ++neoPepMissPos) {
            //
            // Let N be the neo-peptide length, K = [1, ..., N] be
            // the position of the mutation within the neo-peptide,
            // and P be the (unit-offset) position of the mutation
            // in the original protein.  Then, the lower index of
            // the neo-peptide is P - (K - 1).
            //
            int neoPepLower = proteinChangePosition - neoPepMissPos + 1;
            int neoPepUpper = neoPepLower + peptideLength - 1;

            if (neoPepLower < 1 || neoPepUpper > nativeProtein.length())
                continue;

            UnitIndexRange neoPepRange =
                UnitIndexRange.instance(neoPepLower, neoPepUpper);

            NeoPeptide neoPeptide = NeoPeptide.instance(mutatedProtein.fragment(neoPepRange));
            SelfPeptide selfPeptide = SelfPeptide.instance(nativeProtein.fragment(neoPepRange));
            
            Probability neoCleaveProb = mutatedChopEngine.computeCleavageProb(neoPepRange);
            Probability selfCleaveProb = nativeChopEngine.computeCleavageProb(neoPepRange);

            MissCleavageRecord missCleavageRecord =
                MissCleavageRecord.create(tumorBarcode,
                                      hugoSymbol,
                                      proteinChange,
                                      UnitIndex.instance(neoPepMissPos),
                                      neoPepRange,
                                      neoPeptide,
                                      selfPeptide,
                                      neoCleaveProb,
                                      selfCleaveProb);

            missCleavageRecords.add(missCleavageRecord);
        }

        return missCleavageRecords;
    }
}
