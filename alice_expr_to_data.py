#!/usr/bin/env python

import gzip
import sys

from csv import DictReader
from subprocess import Popen, PIPE

from adipo_sight import *
from KGDB import *

fieldnames =["gene",
             "testid",
             "geneid",
             "locus",
             "sample1",
             "sample2",
             "tnf.status",
             "tnf.val1",
             "tnf.val2",
             "tnf.lnFC",
             "tnf.stat",
             "tnf.pval",
             "tnf.qval",
             "tnf.sig",
             "tnf.log2fc",
             "dex.status",
             "dex.val1",
             "dex.val2",
             "dex.lnFC",
             "dex.stat",
             "dex.pval",
             "dex.qval",
             "dex.sig",
             "dex.log2fc",
             "hi.status",
             "hi.val1",
             "hi.val2",
             "hi.lnFC",
             "hi.stat",
             "hi.pval",
             "hi.qval",
             "hi.sig",
             "hi.log2fc",
             "hypo.status",
             "hypo.val1",
             "hypo.val2",
             "hypo.lnFC",
             "hypo.stat",
             "hypo.pval",
             "hypo.qval",
             "hypo.sig",
             "hypo.log2fc",
             "Entrez",
            ]

# locus format = chrX:\d-\d


fn = 'vitro_entrez_all.txt'
f = DictReader(open(fn),delimiter='\t')

kg_sess = get_session('knownGene.db')

# the database is structured as follows:
#
#   - a gene is given a RegionSet record named by the gene
#   - each gene's RegionSet is assigned regions corresponding to the different
#     datatypes, i.e. expression, hypersensitive regions, binding regions,
#     reference regions
#   - each region is assigned a RegionData or SeqData record corresponding to
#     that region's datapoints across conditions


populate_lookup_tables('adipo_sight.db')
adipo_sess = get_session('adipo_sight.db')
found, not_found = 0,0
not_found_list = []

data_types = dict([(r.name,r) for r in adipo_sess.query(DataType).all()])

for data_i, r in enumerate(f) :
    if (data_i % 1000) == 0 :
        print data_i
        adipo_sess.commit()
    gene = r['gene']

    # add gene set
    add_new_or_pass(adipo_sess,[{'name':gene}],RegionSet)
    gene_set = adipo_sess.query(RegionSet).filter(RegionSet.name==gene).first()

    # add region
    chrm, st_en = r['locus'].split(':')
    st, en = map(int, st_en.split('-'))

    region_type = adipo_sess.query(RegionType).filter(RegionType.name=='other').first()
    chrom = adipo_sess.query(Chromosome).filter(Chromosome.name==chrm).first()
    if chrom is None :
        chrom = Chromosome(name=chrm)
        adipo_sess.add(chrom)
        adipo_sess.commit()

    add_d = {'name': '%s diff expr'%r['gene'],
             'region_type_id': region_type.id,
             'chrom_id': chrom.id,
             'start': st,
             'end': en,
             'strand': None,
             'notes': 'differentially expressed region, determined by cufflinks'
            }

    region_recs = add_new_or_pass(adipo_sess,[add_d],Region)
    region_rec = region_recs[0]

    if region_rec is None :
        region_rec = adipo_sess.query(Region) \
                               .filter(Region.name == add_d['name'] and
                                       Region.region_type_id == add_d['region_type_id']).first()
    else :
        region_rec.region_sets.append(gene_set)

    conditions = 'tnf','dex','hi','hypo'

    # val1 is the control condition expression value
    # val2 is the experiment condition expression value
    # log2fc is the log fold change of experiment/condition
    for cond in conditions :
        cond_rec = get_lu_or_add(adipo_sess,cond,Condition)
        for typ, fld in zip(('control expression','experiment expression','log2 expression fold change'),('val1','val2','log2fc')) :
            cond_field = '%s.%s'%(cond,fld)
            data_field_d = {'id':data_i,
                            'region_id':region_rec.id,
                            'data_type_id':data_types[typ].id,
                            'condition_id':cond_rec.id,
                            'value': float(r[cond_field]),
                            'meta1lbl': 'pval',
                            'meta1': r['%s.pval'%cond],
                            'meta2lbl': 'qval',
                            'meta2': r['%s.qval'%cond]
                            }

            adipo_sess.add(RegionData(**data_field_d))
print
