
package pepmhc.junit;

import java.util.List;

import pepmhc.chop.NetChopParser;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopParserTest {
    private static final double TOLERANCE = 0.000001;
    private static final String NETCHOP_FILE = "data/test/netchop.out";

    @Test public void testParse() {
        List<Double> scores = NetChopParser.parse(NETCHOP_FILE);
        assertEquals(123, scores.size());

        assertEquals(0.760600, scores.get(0), TOLERANCE);
        assertEquals(0.483380, scores.get(1), TOLERANCE);
        assertEquals(0.088514, scores.get(2), TOLERANCE);

        assertEquals(0.958277, scores.get(121), TOLERANCE);
        assertEquals(0.177359, scores.get(122), TOLERANCE);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopParserTest");
    }
}
