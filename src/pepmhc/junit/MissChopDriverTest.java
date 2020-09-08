
package pepmhc.junit;

import java.io.File;
import java.util.List;

import jene.peptide.ProteinChange;
import jene.tcga.TumorBarcode;

import pepmhc.miss.MissChopDriver;
import pepmhc.miss.MissChopRecord;
import pepmhc.miss.MissChopTable;

import org.junit.*;
import static org.junit.Assert.*;

public class MissChopDriverTest {
    private static final File outputFile = new File("data/test/miss-chop-test.txt");
    private static final String propFile = "data/test/miss-chop.prop";

    @Test public void testRun() {
        MissChopDriver.run(propFile);

        MissChopTable table = MissChopTable.load(outputFile);
        assertEquals(43, table.count());

        outputFile.delete();
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.MissChopDriverTest");
    }
}
