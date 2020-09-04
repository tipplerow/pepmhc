
package pepmhc.junit;

import java.util.List;

import jene.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.chop.NetChopRunner;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopRunnerTest {
    private static final Peptide PEPTIDE =
        Peptide.instance("MAGRSGDNDEELLKAVRIIKILYKSNPYPEPKGSRQARKN");

    @Test public void testChop() {
        if (!NetChop.isInstalled())
            return;

        List<Double> scores = NetChopRunner.chop(PEPTIDE);

        assertEquals(PEPTIDE.length(), scores.size());
        assertEquals(0.7831, scores.get(3), 0.0001);
        assertEquals(0.0255, scores.get(4), 0.0001);
        assertEquals(0.5506, scores.get(11), 0.0001);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopRunnerTest");
    }
}
