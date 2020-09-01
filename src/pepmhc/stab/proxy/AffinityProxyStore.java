
package pepmhc.stab.proxy;

import jam.sql.SQLStore;

import jene.hla.Allele;

import pepmhc.affy.AffinityMethod;

/**
 * Provides a compute-on-demand service and persistent storage for
 * affinity proxy models.
 */
public final class AffinityProxyStore extends SQLStore<AffinityProxyKey, AffinityProxyModel> {
    private static AffinityProxyStore instance = null;

    private AffinityProxyStore() {
        super(AffinityProxyTable.instance());
    }

    /**
     * Returns the affinity-proxy model for a fixed HLA allele and
     * the default affinity prediction method.
     *
     * @param allele the HLA allele which the model describes.
     *
     * @return the affinity-proxy model for the specified HLA allele
     * and the default affinity prediction method.
     */
    public static AffinityProxyModel get(Allele allele) {
        return get(allele, AffinityProxyPredictor.AFFINITY_METHOD);
    }

    /**
     * Returns the affinity-proxy model for a fixed HLA allele and
     * affinity prediction method.
     *
     * @param allele the HLA allele which the model describes.
     *
     * @param method the affinity prediction method used to build the
     * model.
     *
     * @return the affinity-proxy model for the specified HLA allele
     * and affinity prediction method.
     */
    public static AffinityProxyModel get(Allele allele, AffinityMethod method) {
        return instance().get(AffinityProxyKey.instance(allele, method));
    }

    private static synchronized AffinityProxyStore instance() {
        if (instance == null)
            instance = new AffinityProxyStore();

        return instance;
    }

    @Override protected AffinityProxyModel compute(AffinityProxyKey key) {
        return AffinityProxyBuilder.build(key.getAllele(), key.getMethod());
    }
}
