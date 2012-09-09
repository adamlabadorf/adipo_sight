#!/usr/bin/env python

from collections import defaultdict
from csv import DictReader
from cPickle import dumps
from itertools import count

import sqlalchemy
from numpy import loadtxt, load

from adipo_sight.db import *

conditions = 'dex','hi','tnf','hypoxia'

#conditions = 'fake_p65','dex','hi','tnf','hypoxia','7-2_bg','7-2_nobg'

fields = ("knownGeneID",
         "geneSymbol",
         "chrom",
         "chromStart",
         "chromEnd",
         "name",
         "score",
         "strand",
         "thickStart",
         "thickEnd",
         "itemRgb",
         "blockCount",
         "blockSizes",
         "blockStarts",
         "peak loc",
         "dist from feature",
         "score",
         "map type",
         "map subtype")

if __name__ == '__main__' :

    sess = get_session('adipo_sight.db')

    region_type = get_lu_or_add(sess,'hypersensitive peak',RegionType)
    seq_type = get_lu_or_add(sess,'motif scores',SeqType)

    seq_id_counter = count(1)

    gene_db_d = {}
    region_db_d = {}
    for cond in conditions :
        print 'processing',cond

        # db condition lookup value
        cond_rec = get_lu_or_add(sess,cond,Condition)

        # motif scores are motifs x sequences, load and transpose
        motif_fn = '%s_hypersensitive_peaks_motif_scores.txt'%cond

        #peak_motif_scores = loadtxt('dnase/%s'%motif_fn).T
        motif_fn = '%s_hypersensitive_peaks_motif_scores.npy'%cond
        peak_motif_scores = load('dnase/%s'%motif_fn).T
        print peak_motif_scores.shape

        peak_fn = '%s_hypersensitive_peaks.bed'%cond
        #chr start stop name -10*log10(pval)
        peak_f = DictReader(open('dnase/%s'%peak_fn),delimiter='\t',
                            fieldnames=('chrom','chromStart','chromEnd','name','pval'))

        genes_fn = '%s_hypersensitive_peaks_genes_u3kd2k_tss.txt'%cond
        genes_f = DictReader(open('dnase/%s'%genes_fn),delimiter='\t')
        peak_gene_map = defaultdict(set) 
        print 'loading peak-gene mapping'
        for r in genes_f :
            k = '%s:%s-%s'%(r['chrom'],r['chromStart'],r['chromEnd'])
            peak_gene_map[k].add(r['geneSymbol'].lower())

        print 'iterating peaks'
        i = 0
        for peak, motif_scores in zip(peak_f,peak_motif_scores) :

            if i % 1000 == 0 :
                print '%d'%i
            i += 1

            motif_data = dumps(motif_scores.tolist())

            chrom = get_lu_or_add(sess,peak['chrom'],Chromosome)

            peak_key = '%s:%s-%s'%(peak['chrom'],peak['chromStart'],peak['chromEnd'])

            region = region_db_d.get(peak_key)

            if region is None :
                region = Region(name=peak_key,
                                region_type_id=region_type.id,
                                chrom_id=chrom.id,
                                start=int(peak['chromStart']),
                                end=int(peak['chromEnd']),
                                )

                sess.add(region)
                #sess.commit()
                region_db_d[peak_key] = region

            seq_data = SeqData(region=region,
                               seq_type_id=seq_type.id,
                               condition_id=cond_rec.id,
                               value=motif_data,
                              )
            sess.add(seq_data)

            for sym in peak_gene_map[peak_key] :

                region_set = gene_db_d.get(sym)

                if region_set is None :

                    region_set = RegionSet(name=sym)
                    sess.add(region_set)
                    gene_db_d[sym] = region_set

                region_set.regions.append(region)

                """
                """


                #seq_data = SeqData(id=seq_id_counter.next(),
                """
                try :
                    sess.add(seq_data)
                    sess.commit()
                except sqlalchemy.exc.IntegrityError,e :
                    # this will happen if a peak is mapped to different
                    # isoforms of the same gene in the same condition,
                    # we can safely ignore failures
                    print 'skipping redundant sequence data addition for gene %s'%sym
                    print e
                    sess.rollback()
                    break
                """
        print seq_data.region_id
        sess.commit()
        print seq_data.region_id

