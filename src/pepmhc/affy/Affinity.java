
package pepmhc.affy;

import jam.lang.DomainDouble;
import jam.math.DoubleRange;
import jam.vector.VectorView;

/**
 * Represents a peptide-MHC binding affinity, expressed as an IC50
 * concentration in nanomolar units.
 */
public final class Affinity extends DomainDouble implements Comparable<Affinity> {
    /**
     * Valid range for binding affinities.
     */
    public static final DoubleRange RANGE = DoubleRange.POSITIVE;

    /**
     * Creates a new affinity object.
     *
     * @param value the affinity of the event occurring.
     *
     * @throws RuntimeException unless the affinity is valid.
     */
    public Affinity(double value) {
        super(value, RANGE);
    }

    /**
     * Parses the string representation of an affinity.
     *
     * @param s the string representation of an affinity. 
     *
     * @return a new affinity with the value specfied by the input
     * string.
     *
     * @throws IllegalArgumentException unless the string is a valid
     * affinity.
     */
    public static Affinity parse(String s) {
        return valueOf(Double.parseDouble(s));
    }

    /**
     * Validates an affinity.
     *
     * @param value the affinity to validate.
     *
     * @throws IllegalArgumentException unless the affinity is
     * positive.
     */
    public static void validate(double value) {
        if (!RANGE.contains(value))
            throw new IllegalArgumentException("Invalid affinity.");
    }

    /**
     * Marks a {@code double} value as a {@code Affinity}.
     *
     * @param value the affinity value.
     *
     * @return a {@code Affinity} object having the specified
     * affinity value.
     *
     * @throws IllegalArgumentException if the affinity is
     * negative.
     */
    public static Affinity valueOf(double value) {
        return new Affinity(value);
    }

    /**
     * Converts an array of {@code double} values into affinities.
     *
     * @param values the values to convert.
     *
     * @return an array of corresponding {@code Affinity} objects.
     */
    public static Affinity[] valueOf(double... values) {
        return valueOf(VectorView.wrap(values));
    }

    /**
     * Converts a vector view into an array of affinities.
     *
     * @param values the values to convert.
     *
     * @return an array of corresponding {@code Affinity} objects.
     */
    public static Affinity[] valueOf(VectorView values) {
        Affinity[] result = new Affinity[values.length()];

        for (int index = 0; index < values.length(); index++)
            result[index] = valueOf(values.getDouble(index));

        return result;
    }

    @Override public int compareTo(Affinity that) {
        return compare(this, that);
    }
}
