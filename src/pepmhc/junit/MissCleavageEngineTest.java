
package pepmhc.junit;

import java.util.List;

import jam.app.JamEnv;
import jam.math.UnitIndex;
import jam.math.UnitIndexRange;

import jene.ensembl.EnsemblProteinDb;
import jene.hugo.HugoMaster;
import jene.missense.MissenseGroup;
import jene.missense.MissenseTable;
import jene.peptide.Peptide;
import jene.peptide.ProteinChange;
import jene.tcga.TumorBarcode;

import pepmhc.miss.MissCleavageEngine;
import pepmhc.miss.MissCleavageRecord;

import org.junit.*;
import static org.junit.Assert.*;

public class MissCleavageEngineTest {
    @Test public void testGenerate() {
        String jeneHome = JamEnv.getRequired("JENE_HOME");
        HugoMaster hugoMaster = HugoMaster.load(jeneHome + "/data/test/hugo_master_test.tsv");
        EnsemblProteinDb ensemblDb = EnsemblProteinDb.load(jeneHome + "/data/test/ensembl_test2.fa");

        MissCleavageEngine.initialize(hugoMaster, ensemblDb);

        TumorBarcode barcode1 = TumorBarcode.instance("barcode1");
        TumorBarcode barcode2 = TumorBarcode.instance("barcode2");

        MissenseTable missenseTable = MissenseTable.load(jeneHome + "/data/test/ppe_missense.maf");
        List<MissenseGroup> missenseGroups = missenseTable.group(barcode1);

        assertEquals(1, missenseGroups.size());
        assertEquals(3, missenseGroups.get(0).size());

        List<MissCleavageRecord> cleavageRecords =
            MissCleavageEngine.generate(missenseGroups.get(0), 9);

        assertEquals(27, cleavageRecords.size());
        assertEquals(UnitIndexRange.instance(17, 25), cleavageRecords.get(0).getNeoPeptideNativeRange());
        assertEquals(UnitIndexRange.instance(168, 176), cleavageRecords.get(26).getNeoPeptideNativeRange());

        for (MissCleavageRecord cleavageRecord : cleavageRecords) {
            Peptide neoPeptide = cleavageRecord.getNeoPeptide();
            Peptide selfPeptide = cleavageRecord.getSelfPeptide();
            ProteinChange proteinChange = cleavageRecord.getProteinChange();
            UnitIndex neoPeptideMissensePosition = cleavageRecord.getNeoPeptideMissensePosition();

            assertTrue(neoPeptide.residueAt(neoPeptideMissensePosition).equals(proteinChange.getMutated()));
            assertTrue(selfPeptide.residueAt(neoPeptideMissensePosition).equals(proteinChange.getNative()));
        }
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.MissCleavageEngineTest");
    }
}
