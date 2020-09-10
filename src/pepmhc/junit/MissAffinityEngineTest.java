
package pepmhc.junit;

import java.util.List;

import jene.hla.Allele;
import jene.hla.Genotype;
import jene.peptide.Peptide;
import jene.peptide.ProteinChange;
import jene.tcga.TumorBarcode;
import jene.tcga.TumorGenotypeTable;

import pepmhc.affy.AffinityMethod;
import pepmhc.miss.MissAffinityEngine;
import pepmhc.miss.MissAffinityRecord;
import pepmhc.miss.MissCleavageRecord;
import pepmhc.miss.MissCleavageTable;

import org.junit.*;
import static org.junit.Assert.*;

public class MissAffinityEngineTest {
    @Test public void testGenerate1() {
        Allele A0101 = Allele.instance("A0101");
        Allele A0201 = Allele.instance("A0201");
        Allele B0702 = Allele.instance("B0702");
        Allele C0303 = Allele.instance("C0303");
        Genotype genotype = Genotype.instance(A0101, A0201, B0702, C0303);

        MissCleavageTable cleavageTable =
            MissCleavageTable.load("data/test/miss-cleavage.txt");

        TumorBarcode barcode1 = TumorBarcode.instance("barcode1");
        TumorBarcode barcode2 = TumorBarcode.instance("barcode2");

        List<MissCleavageRecord> cleavageRecords =
            cleavageTable.lookup(barcode1);

        assertEquals(27, cleavageRecords.size());

        List<MissAffinityRecord> affinityRecords =
            MissAffinityEngine.generate(genotype,
                                        AffinityMethod.NET_MHC_PAN,
                                        cleavageRecords);

        assertEquals(4 * 27, affinityRecords.size());

        /*
        for (MissAffinityRecord affinityRecord : affinityRecords)
            System.out.println(affinityRecord);

        for (MissAffinityRecord affinityRecord : affinityRecords) {
            Peptide neoPeptide = affinityRecord.getNeoPeptide();
            Peptide selfPeptide = affinityRecord.getSelfPeptide();
            ProteinChange proteinChange = affinityRecord.getProteinChange();
            UnitIndex neoPeptideMissensePosition = affinityRecord.getNeoPeptideMissensePosition();

            assertTrue(neoPeptide.residueAt(neoPeptideMissensePosition).equals(proteinChange.getMutated()));
            assertTrue(selfPeptide.residueAt(neoPeptideMissensePosition).equals(proteinChange.getNative()));
        }
        */
    }

    @Test public void testGenerate2() {
        MissCleavageTable cleavageTable =
            MissCleavageTable.load("data/test/miss-cleavage.txt");

        TumorGenotypeTable genotypeTable =
            TumorGenotypeTable.load("data/test/tumor-patient.txt",
                                    "data/test/patient-genotype.csv");

        List<MissAffinityRecord> affinityRecords =
            MissAffinityEngine.generate(AffinityMethod.NET_MHC_PAN, cleavageTable, genotypeTable);

        assertEquals(5 * 27 + 3 * 16, affinityRecords.size());
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.MissAffinityEngineTest");
    }
}
