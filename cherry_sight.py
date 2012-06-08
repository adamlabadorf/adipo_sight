import cPickle
import os
import re

import cherrypy
import matplotlib.pyplot as mp
import numpy as np

from jinja2 import Environment, FileSystemLoader, Template

from adipo_sight import *
import log

class AdipoSite:

    def __init__(self) :

        # get a database connection and look for genes
        self.sess = get_session('adipo_sight.db')

        self.template_env = Environment(loader=FileSystemLoader('./'))

    def check_session(fn) :
        def f(self,*args) :
            try :
                self.sess.query(RegionSet).first() # this will fail if the session is invalid
            except :
                self.sess = get_session('adipo_sight.db')
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
            region_set = self.sess.query(RegionSet).filter(RegionSet.name==gene).first()
            if region_set is None :
                log.error('could not find record for gene %s'%gene)
                continue
            else :
                if len(region_set.regions) == 0 :
                    log.error('found record for gene %s, but no region, this is wrong'%gene)
                    continue
                self.region_sets.append(region_set)

        templ_d = {}

        diff_expr_d = self.diff_expr()
        templ_d.update(diff_expr_d)

        motif_d = self.motif_enrich()
        templ_d.update(motif_d)

        top_template = self.template_env.get_template('adipo_site_results.html')
        top_content = top_template.render(**templ_d)

        return top_content

    def diff_expr(self) :

        conditions = [r.name for r in self.sess.query(Condition.name).all()]
        header = ['Gene']+conditions
        data = []

        for region_set in self.region_sets :

            data.append([region_set.name])

            for region in region_set.regions :
                log.debug('Diff Expr Region: %s'%region.name)
                for data_rec in region.region_data :
                    log.debug('Diff Expr Region Data: %s (%s)'%(data_rec.value,data_rec.data_type.name))
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
        self.make_heatmap(data_mat,conditions,pass_genes,heatmap_fn,'log2 fold change')

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

    def motif_enrich(self) :

        d = {}

        # get motif scores
        scores = (self.sess.query(Region,SeqData)
                           .join((RegionSet,Region.region_sets))
                           .join(SeqData)
                           .filter(SeqData.seq_type.has(SeqType.name=='motif scores'))
                           .filter(RegionSet.name.in_(gene_list))
                 ).all()

        # TODO - I was just about to get motif scores out of the db and make
        # a heatmap

        score_mat = []
        for region, seqdata in scores :
            motif_scores = cPickle.loads(seqdata.value)
            score_mat.append(motif_scores)
        log.debug('motif score mat: %d x %d'%(len(score_mat),len(score_mat)))

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

    current_dir = os.path.dirname(os.path.abspath(__file__))
    conf = {'/images': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(current_dir, 'images'),
                        'tools.staticdir.content_types': {'png': 'image/png'}
                       }
           }

    cherrypy.quickstart(AdipoSite(), '/', config=conf)
