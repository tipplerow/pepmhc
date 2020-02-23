
package pepmhc.affy.net;

import java.io.BufferedReader;
import java.io.File;
import java.util.List;

import jam.io.IOUtil;

import pepmhc.affy.AffinityRecord;

/**
 * Parses output written by the {@code netMHC} or {@code netMHCpan}
 * affinity predictor programs.
 */
public final class NetMHCParser extends NetParser {
    private NetMHCParser(BufferedReader reader) {
        super(reader);
    }

    /**
     * Parses an output file written by {@code netMHC} or {@code netMHCpan}.
     *
     * @param file the path to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<AffinityRecord> parse(File file) {
        return parse(IOUtil.openReader(file));
    }

    /**
     * Parses an output file written by {@code netMHC} or {@code netMHCpan}.
     *
     * @param fileName the name to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<AffinityRecord> parse(String fileName) {
        return parse(IOUtil.openReader(fileName));
    }

    /**
     * Parses an output stream written by {@code netMHC} or {@code netMHCpan}.
     *
     * @param reader a reader attached to the output stream.
     *
     * @return a list containing every binding record in the output stream.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<AffinityRecord> parse(BufferedReader reader) {
        NetParser parser = new NetMHCParser(reader);
        return parser.parse();
    }

    @Override public int getPeptideFieldIndex() {
        return 2;
    }

    @Override public int getAffinityFieldIndex() {
        return 12;
    }

    @Override public int getPercentileFieldIndex() {
        return 13;
    }
}
