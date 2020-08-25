
package pepmhc.junit;

import jam.junit.NumericTestBase;

import jean.peptide.Residue;

import pepmhc.tap.TAPMatrix;
import pepmhc.tap.TAPPosition;

import org.junit.*;
import static org.junit.Assert.*;

public class TAPMatrixTest extends NumericTestBase {
    @Test public void testLoad() {
        TAPMatrix matrix = TAPMatrix.consensus();

        assertEquals(20, matrix.nrow());
        assertEquals( 9, matrix.ncol());

        assertEquals(Residue.valueOfCode1('A'), matrix.rowKey(0));
        assertEquals(Residue.valueOfCode1('C'), matrix.rowKey(1));
        assertEquals(Residue.valueOfCode1('D'), matrix.rowKey(2));
        assertEquals(Residue.valueOfCode1('E'), matrix.rowKey(3));
        assertEquals(Residue.valueOfCode1('F'), matrix.rowKey(4));
        assertEquals(Residue.valueOfCode1('G'), matrix.rowKey(5));
        assertEquals(Residue.valueOfCode1('H'), matrix.rowKey(6));
        assertEquals(Residue.valueOfCode1('I'), matrix.rowKey(7));
        assertEquals(Residue.valueOfCode1('K'), matrix.rowKey(8));
        assertEquals(Residue.valueOfCode1('L'), matrix.rowKey(9));
        assertEquals(Residue.valueOfCode1('M'), matrix.rowKey(10));
        assertEquals(Residue.valueOfCode1('N'), matrix.rowKey(11));
        assertEquals(Residue.valueOfCode1('P'), matrix.rowKey(12));
        assertEquals(Residue.valueOfCode1('Q'), matrix.rowKey(13));
        assertEquals(Residue.valueOfCode1('R'), matrix.rowKey(14));
        assertEquals(Residue.valueOfCode1('S'), matrix.rowKey(15));
        assertEquals(Residue.valueOfCode1('T'), matrix.rowKey(16));
        assertEquals(Residue.valueOfCode1('V'), matrix.rowKey(17));
        assertEquals(Residue.valueOfCode1('W'), matrix.rowKey(18));
        assertEquals(Residue.valueOfCode1('Y'), matrix.rowKey(19));

        assertEquals(TAPPosition.NTerm1, matrix.colKey(0));
        assertEquals(TAPPosition.NTerm2, matrix.colKey(1));
        assertEquals(TAPPosition.NTerm3, matrix.colKey(2));
        assertEquals(TAPPosition.Pos4,   matrix.colKey(3));
        assertEquals(TAPPosition.Pos5,   matrix.colKey(4));
        assertEquals(TAPPosition.Pos6,   matrix.colKey(5));
        assertEquals(TAPPosition.Pos7,   matrix.colKey(6));
        assertEquals(TAPPosition.Pos8,   matrix.colKey(7));
        assertEquals(TAPPosition.CTerm,  matrix.colKey(8));

        assertDouble(-1.56, matrix.get(0, 0));
        assertDouble(-0.25, matrix.get(0, 1));
        assertDouble(-0.10, matrix.get(0, 2));
        assertDouble( 0.24, matrix.get(0, 3));
        assertDouble(-0.10, matrix.get(0, 4));
        assertDouble( 0.17, matrix.get(0, 5));
        assertDouble( 0.27, matrix.get(0, 6));
        assertDouble( 0.00, matrix.get(0, 7));
        assertDouble( 0.55, matrix.get(0, 8));

        assertDouble( 0.00, matrix.get( 1, 8));
        assertDouble( 1.83, matrix.get( 2, 8));
        assertDouble( 1.58, matrix.get( 3, 8));
        assertDouble(-2.52, matrix.get( 4, 8));
        assertDouble( 1.41, matrix.get( 5, 8));
        assertDouble( 0.55, matrix.get( 6, 8));
        assertDouble(-0.52, matrix.get( 7, 8));
        assertDouble(-0.45, matrix.get( 8, 8));
        assertDouble(-0.94, matrix.get( 9, 8));
        assertDouble(-0.29, matrix.get(10, 8));
        assertDouble( 1.33, matrix.get(11, 8));
        assertDouble(-0.09, matrix.get(12, 8));
        assertDouble( 0.12, matrix.get(13, 8));
        assertDouble(-1.47, matrix.get(14, 8));
        assertDouble( 2.26, matrix.get(15, 8));
        assertDouble( 0.72, matrix.get(16, 8));
        assertDouble(-0.30, matrix.get(17, 8));
        assertDouble(-0.87, matrix.get(18, 8));
        assertDouble(-2.91, matrix.get(19, 8));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.TAPMatrixTest");
    }
}
