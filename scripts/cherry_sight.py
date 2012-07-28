#!/usr/bin/env python

import cPickle
import os
import random
import re

from collections import defaultdict

import cherrypy
import matplotlib.pyplot as mp
import numpy as np

from jinja2 import Environment, PackageLoader, Template

import adipo_sight.log as log

import adipo_sight.db as db
from adipo_sight.mww import mww

from pkg_resources import resource_filename

class AdipoSite:

    def __init__(self) :

        # get a database connection and look for genes
        self.db_name = db_name = resource_filename('adipo_sight','data/adipo_sight.db')
        self.sess = db.get_session(db_name)

        loader = PackageLoader('adipo_sight','data/tmpl')
        self.template_env = Environment(loader=loader)

    def check_session(fn) :
        def f(self,*args) :
            try :
                self.sess.query(db.RegionSet).first() # this will fail if the session is invalid
            except :
                self.sess = db.get_session(self.db_name)
            return fn(self,*args)
        return f

    @cherrypy.expose
    def index(self,gene_list=None) :

        top = ''
        if gene_list is None :
            top = self.adipo_form()
        else :
            top = self.adipo_results(gene_list)

        bottom = self.log()

        template = self.template_env.get_template('adipo_site.html')
        return template.render(top_content=top,bottom_content=bottom)

    def adipo_form(self) :
        log.info('new search')
        top_template = self.template_env.get_template('adipo_site_form.html')
        top_content = top_template.render()
        return top_content

    @check_session
    def adipo_results(self,gene_list_str=None) :

        gene_list_strip = re.sub(r'\s+',' ',gene_list_str)
        log.info('user submitted genes: %s'%gene_list_strip)

        gene_list = gene_list_strip.split()

        # make sure the genes are unique
        gene_list = list(set(gene_list))

        self.region_sets = []
        for gene in gene_list :
            region_set = self.sess.query(db.RegionSet).filter(db.RegionSet.name==gene).first()
            if region_set is None :
                log.error('could not find record for gene %s'%gene)
                continue
            else :
                if len(region_set.regions) == 0 :
                    log.error('found record for gene %s, but no region, this is wrong'%gene)
                    continue
                self.region_sets.append(region_set)

        templ_d = {}

        #diff_expr_d = self.diff_expr()
        #templ_d.update(diff_expr_d)

        motif_d = self.motif_enrich(gene_list)
        templ_d.update(motif_d)

        top_template = self.template_env.get_template('adipo_site_results.html')
        top_content = top_template.render(**templ_d)

        return top_content

    def diff_expr(self) :

        conditions = [r.name for r in self.sess.query(db.Condition.name).all()]
        header = ['Gene']+conditions
        data = []

        for region_set in self.region_sets :

            data.append([region_set.name])

            for region in region_set.regions :
                log.debug('Diff Expr db.Region: %s'%region.name)
                for data_rec in region.region_data :
                    log.debug('Diff Expr db.Region Data: %s (%s)'%(data_rec.value,data_rec.data_type.name))
                    if data_rec.data_type.name == 'log2 expression fold change' :
                        if float(data_rec.meta2) < 0.05 : # meta2 is qvalue
                            data[-1].append(data_rec.value)
                        else :
                            data[-1].append(0.)


#            # condition data
#            for condition in conditions :
#                cond_data = [reg_data for reg_data in region.region_data if reg_data.condition.name == condition]
#                for r in cond_data :
#                    if r.data_type.name == 'log2 expression fold change' :
#                        data[-1].append(r.value)

        templ_d = {'header':header,
                   'content':data
                  }

        pass_genes = [r[0] for r in data]
        data_mat = np.array([r[1:] for r in data])
        heatmap_fn = 'images/heatmap.png'
        templ_d['heatmap_fn'] = heatmap_fn

        # make the log2 fold change heatmap
        #self.make_heatmap(data_mat,conditions,pass_genes,heatmap_fn,'log2 fold change')

        return templ_d

    def make_heatmap(self,data_mat,conditions,genes,fn,title=None) :

        rev_data_mat = data_mat[::-1,:]

        fig = mp.figure()
        fig.subplots_adjust(left=0.2)

        ax = fig.gca()
        pcolor = ax.pcolor(data_mat)
        fig.colorbar(pcolor)

        ax.set_axes('tight')
        ax.set_xticks([x+.5 for x in range(len(conditions))])
        ax.set_xticklabels(conditions[::-1])
        ax.set_yticks([x+.5 for x in range(len(genes))])
        ax.set_yticklabels(genes[::-1])
        ax.set_title(title)

        fig.savefig(fn)

    def motif_enrich(self,gene_list) :

        d = {}

        # TODO: I need to split out the motif results by condition!!!

        # get motif scores
        scores = (self.sess.query(db.Region,db.SeqData)
                           .join((db.RegionSet,db.Region.region_sets))
                           .join(db.SeqData)
                           .filter(db.SeqData.seq_type.has(db.SeqType.name=='motif scores'))
                           .filter(db.RegionSet.name.in_(gene_list))
                 ).all()

        score_mat = []
        gene_names = []
        for region, seqdata in scores :
            motif_scores = cPickle.loads(seqdata.value)
            score_mat.append(motif_scores)
        log.debug('motif score mat: %d x %d'%(len(score_mat),len(score_mat[0])))
        score_mat = np.array(score_mat)

        all_genes = (self.sess.query(db.Region,db.SeqData)
                     .join((db.RegionSet,db.Region.region_sets))
                     .join(db.SeqData)
                     .filter(db.SeqData.seq_type.has(db.SeqType.name=='motif scores'))
                    ).all()

        # compare motif scores of requested genes to all hypersensitive regions
        # in the dataset
        all_scores = np.zeros((len(all_genes),score_mat.shape[1]))
        for i, (region, seqdata) in enumerate(all_genes) :
            scores = cPickle.loads(seqdata.value)
            all_scores[i] = scores
        log.debug('all_scores: %s'%str(all_scores.shape))
        pvals = mww(score_mat.T,all_scores.T,True)
        log.debug('pvals: %s'%str(pvals.shape))

        motif_name_fn = resource_filename('adipo_sight','data/motif_names.txt')
        motif_names = np.array(open(motif_name_fn).readlines())
        thresh_names, thresh_pvals = motif_names[pvals<0.05], pvals[pvals<0.05]
        thresh_imgs = np.array(['images/motif_logos/%03d_motif.png'%i for i in np.where(pvals<0.05)[0]])
        log.debug('np.where(pvals<0.05): %s'%str([i for i in np.where(pvals<0.05)]))
        d['motifs'] = zip(thresh_names[thresh_pvals.argsort()],
                          thresh_pvals[thresh_pvals.argsort()],
                          thresh_imgs[thresh_pvals.argsort()])

        return d

    def make_boxplot(self,vectors,fn,title=None) :
        fig = mp.figure()
        ax = fig.gca()
        ax.boxplot(vectors,sym='')
        ax.set_title(title)
        fig.savefig(fn)

    def log(self) :
        return log.get_log_content()

if __name__ == "__main__" :

    import adipo_sight
    import pkg_resources
    current_dir = os.path.dirname(os.path.abspath(__file__))
    package_dir = pkg_resources.resource_filename(adipo_sight.__name__,'data')
    conf = {'/images': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(package_dir, 'images'),
                        'tools.staticdir.content_types': {'png': 'image/png'}
                       }
           }

    cherrypy.config.update({ 'server.socket_host': '127.0.0.1', # localhost
                             #'server.socket_host': '18.68.8.98', # for forwarding from web server
                             'server.socket_port': 8888
                           })
    cherrypy.quickstart(AdipoSite(), '/', config=conf)
