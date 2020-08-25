
package pepmhc.junit;

import java.util.List;

import jean.peptide.Peptide;

import pepmhc.chop.NetChop;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopTest {
    private static final Peptide PEPTIDE =
        Peptide.instance("MAGRSGDNDEELLKAVRIIKILYKSNPYPEPKGSRQARKN" +
                         "RRRRWRARQRQIDSISERILSTYLGRSTEPVPLQLPPLER" +
                         "LHLDCREDCGTSGTQQSQGVETGVGRPQISVESPVILGSR" +
                         "TKN");

    private static final List<Peptide> FRAGMENTS =
        List.of(Peptide.instance("SGDNDEELL"),
                Peptide.instance("LKAVRIIKI"),
                Peptide.instance("KAVRIIKIL"),
                Peptide.instance("KILYKSNPY"),
                Peptide.instance("KSNPYPEPK"),
                Peptide.instance("ARKNRRRRW"),
                Peptide.instance("QIDSISERI"),
                Peptide.instance("SERILSTYL"),
                Peptide.instance("YLGRSTEPV"),
                Peptide.instance("GRSTEPVPL"),
                Peptide.instance("STEPVPLQL"),
                Peptide.instance("PLQLPPLER"),
                Peptide.instance("PPLERLHLD"));

    @Test public void testChop() {
        if (NetChop.isInstalled())
            assertEquals(FRAGMENTS, NetChop.chop(PEPTIDE, new int[] { 9 }, 0.5));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopTest");
    }
}
