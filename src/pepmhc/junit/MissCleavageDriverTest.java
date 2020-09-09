
package pepmhc.junit;

import java.io.File;
import java.util.List;

import jene.peptide.ProteinChange;
import jene.tcga.TumorBarcode;

import pepmhc.miss.MissCleavageDriver;
import pepmhc.miss.MissCleavageRecord;
import pepmhc.miss.MissCleavageTable;

import org.junit.*;
import static org.junit.Assert.*;

public class MissCleavageDriverTest {
    private static final File outputFile = new File("data/test/miss-cleavage-test.txt");
    private static final String propFile = "data/test/miss-cleavage.prop";

    @Test public void testRun() {
        MissCleavageDriver.run(propFile);

        MissCleavageTable table = MissCleavageTable.load(outputFile);
        assertEquals(43, table.count());

        outputFile.delete();
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.MissCleavageDriverTest");
    }
}
