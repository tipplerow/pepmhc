
package pepmhc.stab;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.hla.Allele;
import jam.hugo.HugoPeptideTable;
import jam.io.IOUtil;
import jam.io.FileUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;
import jam.util.ListUtil;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

public final class AffinityProxyModel {
    private final Allele allele;
    private final PredictionMethod method;

    private final double intercept;
    private final double coefficient;

    private AffinityProxyModel(Allele allele, PredictionMethod method, double intercept, double coefficient) {
        this.allele = allele;
        this.method = method;

        this.intercept = intercept;
        this.coefficient = coefficient;
    }

    /**
     * System property naming the file that contains a HUGO-peptide
     * table that will provide peptide fragments for model building.
     */
    public static final String PROXY_PEPTIDE_FILE_PROPERTY = "pepmhc.stab.proxyPeptideFile";

    /**
     * System property defining the maximum number of peptides to use
     * in model building.
     */
    public static final String PROXY_SAMPLE_SIZE_PROPERTY = "pepmhc.stab.proxySampleSize";

    /**
     * Default value for the maximum number of peptides to use in
     * model building.
     */
    public static final int PROXY_SAMPLE_SIZE_DEFAULT = 100000;

    /**
     * Prefix for the regression intercept in the parameter file.
     */
    public static final String INTERCEPT_PREFIX = "Intercept:";

    /**
     * Prefix for the regression coefficient in the parameter file.
     */
    public static final String COEFFICIENT_PREFIX = "Coefficient:";

    /**
     * Builds an affinity-proxy model for a given HLA allele and
     * affinity prediction method.
     *
     * @param allele the HLA allele of interest.
     *
     * @param method the affinity prediction method to use.
     *
     * @return an affinity-proxy model with parameters determined for
     * the the specified HLA allele and affinity prediction method.
     */
    public static AffinityProxyModel build(Allele allele, PredictionMethod method) {
        SimpleRegression regression = Builder.build(allele, method);

        assert regression.hasIntercept();
        assert regression.getSlope() < 0.0;

        double intercept = regression.getIntercept();
        double coefficient = regression.getSlope();

        return new AffinityProxyModel(allele, method, intercept, coefficient);
    }

    /**
     * Determines whether a specific affinity-proxy model has been
     * built.
     *
     * @param allele the HLA allele for the model.
     *
     * @param method the affinity prediction method for the model.
     *
     * @return {@code true} iff an affinity-proxy model for the
     * specified allele and method has been built and stored.
     */
    public static boolean exists(Allele allele, PredictionMethod method) {
        File file = file(allele, method);
        return file.exists() && file.canRead();
    }

    /**
     * Returns the parameter file for an affinity-proxy model.
     *
     * @param allele the HLA allele for the model.
     *
     * @param method the affinity prediction method for the model.
     *
     * @return the parameter file for the affinity-proxy model for the
     * specified allele and method.
     */
    public static File file(Allele allele, PredictionMethod method) {
        return new File(path(allele, method));
    }

    /**
     * Returns the affinity-proxy prediction engine for a given HLA
     * allele and affinity prediction method.
     *
     * @param allele the HLA allele of interest.
     *
     * @param method the affinity prediction method.
     *
     * @return the affinity-proxy prediction engine for the specified
     * affinity prediction method.
     */
    public static AffinityProxyModel instance(Allele allele, PredictionMethod method) {
        if (exists(allele, method))
            return load(allele, method);

        AffinityProxyModel model = build(allele, method);
        model.store();

        return model;
    }

    /**
     * Loads the affinity-proxy model for a given HLA allele and
     * affinity prediction method.
     *
     * @param allele the HLA allele of interest.
     *
     * @param method the affinity prediction method to use.
     *
     * @return an affinity-proxy model with parameters loaded from the
     * corresponding parameter file.
     *
     * @throws RuntimeException unless a properly formatted parameter
     * file exists.
     */
    public static AffinityProxyModel load(Allele allele, PredictionMethod method) {
        File file = file(allele, method);
        List<String> lines = IOUtil.readLines(file);

        if (lines.size() != 2)
            throw JamException.runtime("Invalid parameter file: [%s].", file);

        double intercept   = parseIntercept(lines);
        double coefficient = parseCoefficient(lines);

        return new AffinityProxyModel(allele, method, intercept, coefficient);
    }

    private static double parseIntercept(List<String> lines) {
        return parseLine(lines.get(0), INTERCEPT_PREFIX);
    }

    private static double parseCoefficient(List<String> lines) {
        return parseLine(lines.get(1), COEFFICIENT_PREFIX);
    }

    private static double parseLine(String line, String prefix) {
        if (!line.startsWith(prefix))
            throw JamException.runtime("Invalid parameter line: [%s].", line);

        return Double.parseDouble(line.substring(prefix.length()));
    }

    /**
     * Returns the full path name of the parameter file for an
     * affinity-proxy model.
     *
     * @param allele the HLA allele for the model.
     *
     * @param method the affinity prediction method for the model.
     *
     * @return the full path name of the parameter file for the
     * affinity-proxy model for the specified allele and method.
     */
    public static String path(Allele allele, PredictionMethod method) {
        String pepHome  = JamEnv.getRequired("PEPMHC_HOME");
        String baseName = String.format("%s_%s.txt", allele.shortKey(), method);

        return FileUtil.join(pepHome, "data", "proxy", baseName);
    }

    /**
     * Computes the half-life that corresponds to a given affinity
     * according to this proxy model.
     *
     * @param affinity the peptide-MHC binding affinity.
     *
     * @return the half-life corresonding to the specified affinity
     * according to this proxy model.
     *
     * @throws IllegalArgumentException unless the affinity is positive.
     */
    public double halfLife(double affinity) {
        //
        // Recall that the regression model is:
        //
        //     log(halfLife) = intercept + coefficient * log(affinity)
        //
        if (affinity <= 0.0)
            throw new IllegalArgumentException("Negative affinity.");

        return Math.exp(intercept + coefficient * Math.log(affinity));
    }

    /**
     * Stores the parameters of affinity proxy model in the standard
     * location defined by the {@code pepmhc} library.
     */
    public void store() {
        PrintWriter writer = IOUtil.openWriter(file(allele, method));

        writer.println(String.format("%s %9.4f", INTERCEPT_PREFIX, intercept));
        writer.println(String.format("%s %7.4f", COEFFICIENT_PREFIX, coefficient));
        writer.close();
    }

    private static final class Builder {
        private final Allele allele;
        private final PredictionMethod method;

        private double[][] sampleData;
        private List<Peptide> peptideList;
        private SimpleRegression regression;

        private Builder(Allele allele, PredictionMethod method) {
            this.allele = allele;
            this.method = method;
        }

        private static SimpleRegression build(Allele allele, PredictionMethod method) {
            Builder builder = new Builder(allele, method);
            return builder.build();
        }

        private SimpleRegression build() {
            selectSamplePeptides();
            generateSampleData();
            createRegression();
            
            return regression;
        }

        private void selectSamplePeptides() {
            String peptideFile = JamProperties.getRequired(PROXY_PEPTIDE_FILE_PROPERTY);
            int    sampleSize  = JamProperties.getOptionalInt(PROXY_SAMPLE_SIZE_PROPERTY, PROXY_SAMPLE_SIZE_DEFAULT);

            HugoPeptideTable peptideTable = HugoPeptideTable.load(peptideFile);
            peptideList = new ArrayList<Peptide>(peptideTable.viewPeptides());

            if (sampleSize < peptideList.size()) {
                ListUtil.shuffle(peptideList);
                peptideList = peptideList.subList(0, sampleSize);
            }
        }

        private void generateSampleData() {
            List<BindingRecord> affinityRecords = Predictor.predict(method, allele, peptideList);
            List<StabilityRecord> stabilityRecords = NetStab.run(allele, peptideList);

            sampleData = new double[peptideList.size()][];

            for (int k = 0; k < peptideList.size(); ++k) {
                //
                // The linear regression model is: log(H) = intercept + coeff * log(A)
                //
                double affinity = affinityRecords.get(k).getAffinity();
                double halfLife = stabilityRecords.get(k).getHalfLife();

                double logA = Math.log(affinity);
                double logH = Math.log(halfLife);

                sampleData[k] = new double[] { logA, logH };
            }
        }

        private void createRegression() {
            regression = new SimpleRegression();
            regression.addData(sampleData);
        }
    }
}
