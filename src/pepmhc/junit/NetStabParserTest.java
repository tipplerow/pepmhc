
package pepmhc.junit;

import java.util.List;

import jene.peptide.Peptide;

import pepmhc.stab.net.NetStabParser;
import pepmhc.stab.StabilityRecord;

import org.junit.*;
import static org.junit.Assert.*;

public class NetStabParserTest {
    private static final String FILE_NAME = "data/test/netMHCstabpan.out";

    @Test public void testParse() {
        List<StabilityRecord> records = NetStabParser.parse(FILE_NAME);

        assertEquals(10, records.size());

        assertEquals(Peptide.instance("AAAWYLWEV"), records.get(0).getPeptide());
        assertEquals(Peptide.instance("AEFGPWQTV"), records.get(9).getPeptide());
        
        assertEquals(7.61, records.get(0).getHalfLife(), 0.01);
        assertEquals(0.50, records.get(9).getHalfLife(), 0.01);

        assertEquals(0.30, records.get(0).getPercentile(), 0.01);
        assertEquals(9.00, records.get(9).getPercentile(), 0.01);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetStabParserTest");
    }
}
