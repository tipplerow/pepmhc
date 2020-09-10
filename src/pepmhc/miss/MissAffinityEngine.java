
package pepmhc.miss;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamLogger;
import jam.math.Percentile;
import jam.util.ListUtil;
import jam.util.StreamUtil;

import jene.hla.Allele;
import jene.hla.Genotype;
import jene.peptide.Peptide;
import jene.tcga.TumorBarcode;
import jene.tcga.TumorGenotypeTable;

import pepmhc.affy.Affinity;
import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityPredictor;
import pepmhc.bind.BindRecord;
import pepmhc.bind.BindRecordMap;

/**
 * Computes the MHC binding affinity for neo-peptide and self-peptide
 * pairs contained in cleavage records.
 */
public final class MissAffinityEngine {
    private final Allele allele;
    private final AffinityMethod affinityMethod;
    private final Collection<MissCleavageRecord> cleavageRecords;

    private BindRecordMap bindingMap;
    private Collection<Peptide> peptides;
    private List<MissAffinityRecord> affinityRecords;

    private MissAffinityEngine(Allele allele,
                               AffinityMethod affinityMethod,
                               Collection<MissCleavageRecord> cleavageRecords) {
        this.allele = allele;
        this.affinityMethod = affinityMethod;
        this.cleavageRecords = cleavageRecords;
    }

    /**
     * Computes MHC binding affinity for neo-peptide and self-peptide
     * pairs contained in cleavage records.
     *
     * @param allele the HLA allele.
     *
     * @param affinityMethod the enumerated binding affinity method.
     *
     * @param cleavageRecords the cleavage records that contain the
     * target peptides.
     *
     * @return a list of affinity records for the specified cleavage
     * records and allele.
     */
    public static List<MissAffinityRecord> generate(Allele allele,
                                                    AffinityMethod affinityMethod,
                                                    Collection<MissCleavageRecord> cleavageRecords) {
        MissAffinityEngine engine =
            new MissAffinityEngine(allele, affinityMethod, cleavageRecords);

        try {
            return engine.generate();
        }
        catch (RuntimeException ex) {
            JamLogger.warn(ex);
            return List.of();
        }
    }

    /**
     * Computes MHC binding affinity for neo-peptide and self-peptide
     * pairs contained in cleavage records.
     *
     * @param genotype the HLA genotype.
     *
     * @param affinityMethod the enumerated binding affinity method.
     *
     * @param cleavageRecords the cleavage records that contain the
     * target peptides.
     *
     * @return a list of affinity records for the specified cleavage
     * records and genotype.
     */
    public static List<MissAffinityRecord> generate(Genotype genotype,
                                                    AffinityMethod affinityMethod,
                                                    Collection<MissCleavageRecord> cleavageRecords) {
        List<MissAffinityRecord> affinityRecords =
            new ArrayList<MissAffinityRecord>();

        for (Allele allele : genotype.viewUniqueAlleles())
            affinityRecords.addAll(generate(allele, affinityMethod, cleavageRecords));

        return affinityRecords;
    }

    /**
     * Computes MHC binding affinity for neo-peptide and self-peptide
     * pairs contained in cleavage records.
     *
     * @param affinityMethod the enumerated binding affinity method.
     *
     * @param cleavageTable the cleavage records that contain the
     * target peptides.
     *
     * @param genotypeTable a table containing the patient genotypes
     * for each tumor in the cleavage records.
     *
     * @return a list of affinity records for the specified cleavage
     * records and tumor genotypes.
     */
    public static List<MissAffinityRecord> generate(AffinityMethod affinityMethod,
                                                    MissCleavageTable cleavageTable,
                                                    TumorGenotypeTable genotypeTable) {
        Collection<TumorBarcode> tumorBarcodes =
            genotypeTable.viewBarcodes();

        List<List<MissAffinityRecord>> tumorRecords =
            StreamUtil.applyParallel(tumorBarcodes,
                                     tumorBarcode -> generate(tumorBarcode,
                                                              affinityMethod,
                                                              cleavageTable,
                                                              genotypeTable));
        List<MissAffinityRecord> affinityRecords =
            ListUtil.cat(tumorRecords);

        JamLogger.info("Sorting affinity records...");
        affinityRecords.sort(MissAffinityRecord.COMPARATOR);

        return affinityRecords;
    }

    private static List<MissAffinityRecord> generate(TumorBarcode tumorBarcode,
                                                     AffinityMethod affinityMethod,
                                                     MissCleavageTable cleavageTable,
                                                     TumorGenotypeTable genotypeTable) {
        List<MissCleavageRecord> cleavageRecords =
            cleavageTable.lookup(tumorBarcode);

        if (cleavageRecords.isEmpty())
            return List.of();

        Genotype genotype =
            genotypeTable.lookup(tumorBarcode);

        if (genotype != null)
            return generate(genotype, affinityMethod, cleavageRecords);

        JamLogger.warn("Missing genotype for [%s].", tumorBarcode);
        return List.of();
    }


    private List<MissAffinityRecord> generate() {
        JamLogger.info("Generating [%d] affinity records: [%s, %s]...",
                       cleavageRecords.size(), allele, affinityMethod);

        collectPeptides();
        computeAffinity();
        createRecords();
        
        return affinityRecords;
    }

    private void collectPeptides() {
        peptides = MissCleavageRecord.extractPeptides(cleavageRecords);
    }

    private void computeAffinity() {
        bindingMap = affinityMethod.getPredictor().map(allele, peptides);
    }

    private void createRecords() {
        affinityRecords = new ArrayList<MissAffinityRecord>();

        for (MissCleavageRecord cleavageRecord : cleavageRecords)
            affinityRecords.add(createRecord(cleavageRecord));
    }

    private MissAffinityRecord createRecord(MissCleavageRecord cleavageRec) {
        Peptide neoPeptide = cleavageRec.getNeoPeptide();
        Peptide selfPeptide = cleavageRec.getSelfPeptide();

        BindRecord neoRecord = bindingMap.require(neoPeptide);
        BindRecord selfRecord = bindingMap.require(selfPeptide);

        Affinity neoAffinity = neoRecord.getAffinity();
        Affinity selfAffinity = selfRecord.getAffinity();

        Percentile neoAffinityRank = neoRecord.getPercentile();
        Percentile selfAffinityRank = selfRecord.getPercentile();
        
        return MissAffinityRecord.create(cleavageRec,
                                         allele,
                                         neoAffinity,
                                         selfAffinity,
                                         neoAffinityRank,
                                         selfAffinityRank);
    }
}
