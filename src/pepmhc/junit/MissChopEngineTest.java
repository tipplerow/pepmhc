
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

import pepmhc.miss.MissChopEngine;
import pepmhc.miss.MissChopRecord;

import org.junit.*;
import static org.junit.Assert.*;

public class MissChopEngineTest {
    @Test public void testGenerate() {
        String jeneHome = JamEnv.getRequired("JENE_HOME");
        HugoMaster hugoMaster = HugoMaster.load(jeneHome + "/data/test/hugo_master_test.tsv");
        EnsemblProteinDb ensemblDb = EnsemblProteinDb.load(jeneHome + "/data/test/ensembl_test2.fa");

        MissChopEngine.initialize(hugoMaster, ensemblDb);

        TumorBarcode barcode1 = TumorBarcode.instance("barcode1");
        TumorBarcode barcode2 = TumorBarcode.instance("barcode2");

        MissenseTable missenseTable = MissenseTable.load(jeneHome + "/data/test/ppe_missense.maf");
        List<MissenseGroup> missenseGroups = missenseTable.group(barcode1);

        assertEquals(1, missenseGroups.size());
        assertEquals(3, missenseGroups.get(0).size());

        List<MissChopRecord> chopRecords =
            MissChopEngine.generate(missenseGroups.get(0), 9);

        assertEquals(27, chopRecords.size());
        assertEquals(UnitIndexRange.instance(17, 25), chopRecords.get(0).getNeoPeptideNativeRange());
        assertEquals(UnitIndexRange.instance(168, 176), chopRecords.get(26).getNeoPeptideNativeRange());

        for (MissChopRecord chopRecord : chopRecords) {
            Peptide neoPeptide = chopRecord.getNeoPeptide();
            Peptide selfPeptide = chopRecord.getSelfPeptide();
            ProteinChange proteinChange = chopRecord.getProteinChange();
            UnitIndex neoPeptideMissensePosition = chopRecord.getNeoPeptideMissensePosition();

            assertTrue(neoPeptide.residueAt(neoPeptideMissensePosition).equals(proteinChange.getMutated()));
            assertTrue(selfPeptide.residueAt(neoPeptideMissensePosition).equals(proteinChange.getNative()));
        }

        /*
        assertEquals(UnitIndexRange.instance(9, 17), chopRecords.get(0).getPeptideRange());
        assertEquals(UnitIndexRange.instance(20, 28), chopRecords.get(11).getPeptideRange());
        assertEquals(UnitIndexRange.instance(168, 176), chopRecords.get(12).getPeptideRange());
        assertEquals(UnitIndexRange.instance(176, 184), chopRecords.get(20).getPeptideRange());

        missenseGroups = missenseTable.group(barcode2);
        chopRecords = MissChopEngine.generate(missenseGroups.get(0), 9);

        assertEquals(1, missenseGroups.size());
        assertEquals(16, chopRecords.size());

        assertEquals(UnitIndexRange.instance(1, 9), chopRecords.get(0).getPeptideRange());
        assertEquals(UnitIndexRange.instance(3, 11), chopRecords.get(2).getPeptideRange());
        assertEquals(UnitIndexRange.instance(168, 176), chopRecords.get(3).getPeptideRange());
        assertEquals(UnitIndexRange.instance(181, 189), chopRecords.get(15).getPeptideRange());

        chopRecords = MissChopEngine.generate(missenseTable, 9);

        assertEquals(37, chopRecords.size());

        assertEquals(barcode1, chopRecords.get(0).getTumorBarcode());
        assertEquals(UnitIndexRange.instance(9, 17), chopRecords.get(0).getPeptideRange());

        assertEquals(barcode2, chopRecords.get(36).getTumorBarcode());
        assertEquals(UnitIndexRange.instance(181, 189), chopRecords.get(36).getPeptideRange());
        */
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.MissChopEngineTest");
    }
}
