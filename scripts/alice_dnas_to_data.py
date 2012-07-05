#!/usr/bin/env python

from csv import DictReader
from cPickle import dumps
from itertools import count

import sqlalchemy
from numpy import loadtxt

from adipo_sight.db import *

conditions = 'dex','hi','tnf','hypoxia'

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

    for cond in conditions :

        # db condition lookup value
        cond_rec = get_lu_or_add(sess,cond,Condition)

        # motif scores are motifs x sequences, load and transpose
        motif_fn = '%s_hypersensitive_peaks_motif_scores.txt'%cond
        peak_motif_scores = loadtxt('dnase/%s'%motif_fn).T

        genes_fn = '%s_hypersensitive_peaks_genes.txt'%cond
        genes_f = DictReader(open('dnase/%s'%genes_fn),delimiter='\t')
        gene_db_d = {}
        for peak, motif_scores in zip(genes_f,peak_motif_scores) :

            sym = peak['geneSymbol']

            if sym not in gene_db_d :

                region_set_recs = sess.query(RegionSet).filter(RegionSet.name == sym).all()

                if len(region_set_recs) == 0 :
                    # this shouldn't happen, but no RegionSet found
                    continue

                elif len(region_set_recs) > 1 :
                    print 'this should never happen'

                gene_db_d[sym] = region_set_recs[0]

            region_set = gene_db_d[sym]

            chrom = get_lu_or_add(sess,peak['chrom'],Chromosome)

            region = Region(name='%s DNase %s:%s-%s'%(sym,peak['chrom'],
                                                      peak['chromStart'],
                                                      peak['chromEnd']
                                                     ),
                            region_type_id=region_type.id,
                            chrom_id=chrom.id,
                            start=int(peak['chromStart']),
                            end=int(peak['chromEnd']),
                            )

            region_set.regions.append(region)

            try :
                sess.add(region)
                sess.commit()
            except sqlalchemy.exc.IntegrityError :
                # this will happen if a peak is mapped to different
                # isoforms of the same gene, we can safely ignore failures
                print 'skipping redundant region addition for gene %s'%sym
                sess.rollback()
                region = sess.query(Region).filter(Region.name==region.name,Region.region_type_id==region_type.id).first()

            motif_data = dumps(motif_scores.tolist())

            seq_data = SeqData(id=seq_id_counter.next(),
                               region_id=region.id,
                               seq_type_id=seq_type.id,
                               condition_id=cond_rec.id,
                               value=motif_data,
                               meta1lbl='dist from feature',
                               meta1=peak['dist from feature']
                              )
            try :
                sess.add(seq_data)
                sess.commit()
            except sqlalchemy.exc.IntegrityError :
                # this will happen if a peak is mapped to different
                # isoforms of the same gene in the same condition,
                # we can safely ignore failures
                print 'skipping redundant sequence data addition for gene %s'%sym
                sess.rollback()

