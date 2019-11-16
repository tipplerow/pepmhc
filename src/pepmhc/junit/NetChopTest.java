
package pepmhc.junit;

import java.util.List;

import jam.peptide.Peptide;
import pepmhc.chop.NetChop;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopTest {
    private static final Peptide PEPTIDE =
        Peptide.parse("MAGRSGDNDEELLKAVRIIKILYKSNPYPEPKGSRQARKN" +
                      "RRRRWRARQRQIDSISERILSTYLGRSTEPVPLQLPPLER" +
                      "LHLDCREDCGTSGTQQSQGVETGVGRPQISVESPVILGSR" +
                      "TKN");

    private static final List<Peptide> FRAGMENTS =
        List.of(Peptide.parse("SGDNDEELL"),
                Peptide.parse("LKAVRIIKI"),
                Peptide.parse("KAVRIIKIL"),
                Peptide.parse("KILYKSNPY"),
                Peptide.parse("KSNPYPEPK"),
                Peptide.parse("ARKNRRRRW"),
                Peptide.parse("QIDSISERI"),
                Peptide.parse("SERILSTYL"),
                Peptide.parse("YLGRSTEPV"),
                Peptide.parse("GRSTEPVPL"),
                Peptide.parse("STEPVPLQL"),
                Peptide.parse("PLQLPPLER"),
                Peptide.parse("PPLERLHLD"));

    @Test public void testChop() {
        if (NetChop.isInstalled())
            assertEquals(FRAGMENTS, NetChop.chop(PEPTIDE, new int[] { 9 }, 0.5));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopTest");
    }
}
