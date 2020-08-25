
package pepmhc.stab.proxy;

import java.util.regex.Pattern;

import jam.lang.JamException;
import jam.util.RegexUtil;

import jean.hla.Allele;

import pepmhc.affy.AffinityMethod;

/**
 * Uniquely identifies an affinity proxy model.
 */
public final class AffinityProxyKey {
    private final Allele allele;
    private final AffinityMethod method;
    private final String string;

    private static final Pattern PATTERN   = RegexUtil.PIPE;
    private static final String  SEPARATOR = "|";

    private AffinityProxyKey(Allele allele, AffinityMethod method) {
        this.allele = allele;
        this.method = method;
        this.string = format(allele, method);
    }

    private static String format(Allele allele, AffinityMethod method) {
        return allele.shortKey() + SEPARATOR + method.name();
    }

    /**
     * Returns the key for a given allele and affinity method.
     *
     * @param allele the allele for the key.
     *
     * @param method the affinity prediction method for the key.
     *
     * @return the key for the specified allele and affinity method.
     */
    public static AffinityProxyKey instance(Allele allele, AffinityMethod method) {
        return new AffinityProxyKey(allele, method);
    }

    /**
     * Returns the key encoded by a string.
     *
     * @param string the string representation to parse.
     *
     * @return the key encoded by the specified string.
     */
    public static AffinityProxyKey parse(String string) {
        String[] fields = RegexUtil.split(PATTERN, string, 2);

        Allele allele = Allele.instance(fields[0]);
        AffinityMethod method = AffinityMethod.valueOf(fields[1]);

        return instance(allele, method);
    }

    /**
     * Formats this key for storage in the parameter database.
     *
     * @return the string representation of this key stored in
     * the parameter database.
     */
    public String format() {
        return string;
    }

    /**
     * Returns the allele for this key.
     *
     * @return the allele for this key.
     */
    public Allele getAllele() {
        return allele;
    }

    /**
     * Returns the affinity prediction method for this key.
     *
     * @return the affinity prediction method for this key.
     */
    public AffinityMethod getMethod() {
        return method;
    }

    @Override public boolean equals(Object obj) {
        return (obj instanceof AffinityProxyKey) && equalsKey((AffinityProxyKey) obj);
    }

    private boolean equalsKey(AffinityProxyKey that) {
        return this.string.equals(that.string);
    }

    @Override public int hashCode() {
        return string.hashCode();
    }

    @Override public String toString() {
        return string;
    }
}
