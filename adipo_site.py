#!/usr/bin/python

import cgi
import cgitb
cgitb.enable()

import Cookie
import numpy as np
import os
import re
import sha
import sys
import time

import matplotlib.pyplot as mp

from jinja2 import Environment, FileSystemLoader, Template

from adipo_sight import *
import log

# set up logging
print "Content-Type: text/html"
print

env = Environment(loader=FileSystemLoader('./'))
template = env.get_template('adipo_site.html')

form = cgi.FieldStorage()
gene_list_str = form.getvalue("gene_list")

# this happens when there is a query
if gene_list_str is not None :

    # parse gene list, collapse all contiguous whitespace to one space
    gene_list_strip = re.sub(r'\s+',' ',gene_list_str)
    log.info('user submitted genes: %s'%gene_list_strip)

    gene_list = gene_list_strip.split()

    # get a database connection and look for genes
    sess = get_session('adipo_sight.db')

    conditions = [r.name for r in sess.query(Condition.name).all()]
    header = ['Gene']+conditions
    data = []
    for gene in gene_list :
        region_set = sess.query(RegionSet).filter(RegionSet.name == gene).first()
        if region_set is None :
            log.error('could not find record for gene %s, this is wrong'%gene)
            continue
        else :
            region = region_set.regions[0] # only one region per gene right now

        data.append([gene])

        # condition data
        for condition in conditions :
            cond_data = [reg_data for reg_data in region.region_data if reg_data.condition.name == condition]
            for r in cond_data :
                if r.data_type.name == 'log2 expression fold change' :
                    data[-1].append(r.value)

    pass_genes = [r[0] for r in data]
    data_mat = np.array([r[1:] for r in data])
    heatmap_fn = 'images/heatmap.png'

    mp.subplots_adjust(left=0.2)
    mp.pcolor(data_mat)
    mp.colorbar()
    mp.xticks([x+.5 for x in range(len(conditions))],conditions)
    mp.yticks([x+.5 for x in range(len(pass_genes))],pass_genes)
    mp.axis('tight')
    mp.savefig(heatmap_fn)

    top_template = env.get_template('adipo_site_results.html')
    top_content = top_template.render(header=header,content=data,heatmap_fn=heatmap_fn)
    sess.close()

# this is a new search
else :
    log.info('new search')
    top_template = env.get_template('adipo_site_form.html')
    top_content = top_template.render()


bottom_content = log.get_log_content() 
print template.render(top_content=top_content,bottom_content=bottom_content)
