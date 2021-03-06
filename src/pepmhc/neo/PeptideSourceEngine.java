
package pepmhc.neo;

import java.util.List;
import java.util.stream.Collectors;

import jam.app.JamLogger;
import jam.lang.JamException;

import jene.hugo.HugoPeptideList;
import jene.hugo.HugoPeptideTable;
import jene.hugo.HugoSymbol;
import jene.maf.MAFFastaList;
import jene.maf.MAFFastaRecord;
import jene.peptide.Peptide;
import jene.tcga.TumorBarcode;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import pepmhc.agpro.AntigenProcessor;

final class PeptideSourceEngine {
    private final TumorBarcode barcode;
    private final MAFFastaList fastaList;
    private final HugoPeptideTable selfPepReference;
    private final String agproPropFile;

    private final Multimap<HugoSymbol, Peptide> neoPeptideMap = HashMultimap.create();
    private final Multimap<HugoSymbol, Peptide> selfPeptideMap = HashMultimap.create();

    private int processed = 0;

    private PeptideSourceEngine(TumorBarcode barcode,
                                MAFFastaList fastaList,
                                HugoPeptideTable selfPepReference,
                                String agproPropFile) {
        this.barcode = barcode;
        this.fastaList = fastaList;
        this.agproPropFile = agproPropFile;
        this.selfPepReference = selfPepReference;
    }

    static PeptideSourceView process(TumorBarcode barcode,
                                     MAFFastaList fastaList,
                                     HugoPeptideTable selfPepReference,
                                     String agproPropFile) {
        PeptideSourceEngine engine =
            new PeptideSourceEngine(barcode, fastaList, selfPepReference, agproPropFile);

        return engine.process();
    }

    private PeptideSourceView process() {
        processRecords();
        
        HugoPeptideTable neoPeptideTable = HugoPeptideTable.create(neoPeptideMap);
        HugoPeptideTable selfPeptideTable = HugoPeptideTable.create(selfPeptideMap);

        return PeptideSourceView.create(neoPeptideTable, selfPeptideTable, selfPepReference);
    }

    private void processRecords() {
        List<HugoPeptideList> fragmentLists =
            fastaList.parallelStream().map(record -> processRecord(record)).collect(Collectors.toList());

        for (HugoPeptideList fragmentList : fragmentLists)
            processFragmentList(fragmentList);
    }

    private HugoPeptideList processRecord(MAFFastaRecord record) {
        ++processed;

        HugoSymbol symbol = record.getHugoSymbol();
        JamLogger.info("Antigen processing [%s:%s] (%d of %d)...",
                       barcode.getKey(), symbol.getKey(),
                       processed, fastaList.size());
        
        AntigenProcessor antigenProcessor =
            agproPropFile != null ?
            AntigenProcessor.resolve(agproPropFile) : 
            AntigenProcessor.defaultProcessor();

        Peptide peptide = record.getPeptide();
        List<Peptide> fragments = antigenProcessor.process(peptide);

        return HugoPeptideList.wrap(symbol, fragments);
    }

    private void processFragmentList(HugoPeptideList fragmentList) {
        HugoSymbol symbol = fragmentList.getSymbol();
        JamLogger.info("Processing fragment list [%s:%s]...", barcode.getKey(), symbol.getKey());

        for (Peptide fragment : fragmentList) {
            if (isSelfPeptide(fragment))
                selfPeptideMap.put(symbol, fragment);
            else
                neoPeptideMap.put(symbol, fragment);
        }
    }

    private boolean isSelfPeptide(Peptide peptide) {
        return selfPepReference.contains(peptide);
    }
}
