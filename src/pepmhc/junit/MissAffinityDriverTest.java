
package pepmhc.junit;

import java.io.File;

import pepmhc.miss.MissAffinityDriver;
import pepmhc.miss.MissAffinityRecord;
import pepmhc.miss.MissAffinityTable;

import org.junit.*;
import static org.junit.Assert.*;

public class MissAffinityDriverTest {
    private static final File outputFile = new File("data/test/miss-affinity-test.txt");
    private static final String propFile = "data/test/miss-affinity.prop";

    @Test public void testRun() {
        MissAffinityDriver.run(propFile);

        MissAffinityTable table = MissAffinityTable.load(outputFile);
        assertEquals(183, table.count());

        outputFile.delete();
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.MissAffinityDriverTest");
    }
}
