
package pepmhc.engine.net;

import java.io.BufferedReader;
import java.io.File;
import java.util.List;

import jam.io.IOUtil;

import pepmhc.binder.BindingRecord;

/**
 * Parses output written by the {@code netMHCstabpan} stability
 * predictor program.
 */
public final class NetStabParser extends NetParser {
    private NetStabParser(BufferedReader reader) {
        super(reader);
    }

    /**
     * Parses an output file written by {@code netMHCstabpan}.
     *
     * @param file the path to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<BindingRecord> parse(File file) {
        return parse(IOUtil.openReader(file));
    }

    /**
     * Parses an output file written by {@code netMHCstabpan}.
     *
     * @param fileName the name to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<BindingRecord> parse(String fileName) {
        return parse(IOUtil.openReader(fileName));
    }

    /**
     * Parses an output stream written by {@code netMHCstabpan}.
     *
     * @param reader a reader attached to the output stream.
     *
     * @return a list containing every binding record in the output
     * stream.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<BindingRecord> parse(BufferedReader reader) {
        NetParser parser = new NetStabParser(reader);
        return parser.parse();
    }

    @Override public int getPeptideFieldIndex() {
        return 2;
    }

    @Override public int getAffinityFieldIndex() {
        return 5;
    }

    @Override public int getPercentileFieldIndex() {
        return 6;
    }
}
