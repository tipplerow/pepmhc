
package pepmhc.tap;

/**
 * Distinguishes the residue positions that contribute to TAP binding
 * predictions.
 *
 * <p>See the consensus scoring matrix in Table 1 of Peters et al.,
 * <em>J. Immunol.</em> <b>171</b>, 1741-1749 (2003).
 */
public enum TAPPosition {
    NTerm1,
    NTerm2,
    NTerm3,
    Pos4,
    Pos5,
    Pos6,
    Pos7,
    Pos8,
    CTerm;
}
