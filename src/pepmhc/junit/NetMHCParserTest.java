
package pepmhc.junit;

import java.util.List;

import jene.peptide.Peptide;

import pepmhc.affy.AffinityRecord;
import pepmhc.affy.net.NetMHCParser;

import org.junit.*;
import static org.junit.Assert.*;

public class NetMHCParserTest {
    private static final String NET_MHC_FILE = "data/test/netMHC.out";
    private static final String NET_MHC_PAN_FILE = "data/test/netMHCpan.out";

    @Test public void testNetMHC() {
        List<AffinityRecord> records = NetMHCParser.parse(NET_MHC_FILE);

        assertEquals(10, records.size());

        assertEquals(Peptide.instance("AAAWYLWEV"), records.get(0).getPeptide());
        assertEquals(Peptide.instance("AEFGPWQTV"), records.get(9).getPeptide());
        
        assertEquals(   10.65, records.get(0).getAffinity().doubleValue(), 0.01);
        assertEquals(12361.73, records.get(9).getAffinity().doubleValue(), 0.01);

        assertEquals( 0.12, records.get(0).getPercentile().doubleValue(), 0.01);
        assertEquals(18.00, records.get(9).getPercentile().doubleValue(), 0.01);
    }

    @Test public void testNetMHCPan() {
        List<AffinityRecord> records = NetMHCParser.parse(NET_MHC_PAN_FILE);

        assertEquals(10, records.size());

        assertEquals(Peptide.instance("AAAWYLWEV"), records.get(0).getPeptide());
        assertEquals(Peptide.instance("AEFGPWQTV"), records.get(9).getPeptide());
        
        assertEquals(   7.1, records.get(0).getAffinity().doubleValue(), 0.1);
        assertEquals(8706.7, records.get(9).getAffinity().doubleValue(), 0.1);

        assertEquals( 0.0687, records.get(0).getPercentile().doubleValue(), 0.0001);
        assertEquals(13.2114, records.get(9).getPercentile().doubleValue(), 0.0001);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetMHCParserTest");
    }
}
