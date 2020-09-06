
package pepmhc.junit;

import java.util.List;

import jam.math.Probability;

import pepmhc.chop.NetChopParser;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopParserTest {
    private static final double TOLERANCE = 0.000001;
    private static final String NETCHOP_FILE = "data/test/netchop.out";

    @Test public void testParse() {
        List<Probability> scores = NetChopParser.parse(NETCHOP_FILE);
        assertEquals(123, scores.size());

        assertTrue(scores.get(0).equals(0.760600, TOLERANCE));
        assertTrue(scores.get(1).equals(0.483380, TOLERANCE));
        assertTrue(scores.get(2).equals(0.088514, TOLERANCE));

        assertTrue(scores.get(121).equals(0.958277, TOLERANCE));
        assertTrue(scores.get(122).equals(0.177359, TOLERANCE));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopParserTest");
    }
}
