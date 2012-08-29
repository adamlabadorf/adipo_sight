#!/usr/bin/env python

import cPickle
import logging
import os
import random
import re
import sys
import threading

from collections import defaultdict
from itertools import cycle

import cherrypy
import matplotlib.pyplot as mp
import numpy as np

from jinja2 import Environment, PackageLoader, Template

import adipo_sight.log as log

import adipo_sight.db as db
from adipo_sight.mww import mww, mww_multiprocess

from pkg_resources import resource_filename

# some jinja2 custom filters
def dictofdictsort(d,subkey) :
    itms = d.items()
    if subkey is not None :
        itms.sort(key=lambda x: x[1].get(subkey,sys.maxint))
    return itms

class AdipoThread(object) :

    def __init__(self) :
        self.thread = None

    def run_thread(self,target,*args,**kwargs) :

        if self.is_running() :
            raise Exception('Only one thread may be run at a time')

        self.thread = threading.Thread(target=target,args=args,kwargs=kwargs)
        self.thread.start()

    def is_running(self) :
        return self.thread is not None and self.thread.is_alive()

    def is_done(self) :
        is_done = self.thread is not None and not self.thread.is_alive()
        if is_done :
            self.thread = None

class AdipoSite:

    def __init__(self) :

        # get a database connection and look for genes
        self.db_name = db_name = resource_filename('adipo_sight','data/adipo_sight.db')
        cherrypy.log('loading db session from %s'%self.db_name)
        self.db_session = db.get_session(db_name)

        loader = PackageLoader('adipo_sight','data/tmpl')
        self.template_env = Environment(loader=loader)
        self.template_env.filters['dictofdictsort'] = dictofdictsort
        self.template_env.globals['zip'] = zip

        self.args = {}

        self.threads = defaultdict(AdipoThread)

    def check_db_session(fn) :
        def f(self,*args) :
            try :
                self.db_session.query(db.RegionSet).first() # this will fail if the session is invalid
            except :
                self.db_session = db.get_session(self.db_name)
            return fn(self,*args)
        return f

    def parse_kwargs(self,kwargs) :

        for k,v in kwargs.items() :
            if v == '' :
                kwargs[k] = None

        self.args = {'diff_exp_pval':0.05,
                     'diff_exp_hi':1.2,
                     'diff_exp_low':-1.2,
                     'diff_hyp_pval':0.05}
        self.args.update(kwargs)

    @cherrypy.expose
    def index(self,**kwargs) :

        cherrypy.log('loading index')

        self.parse_kwargs(kwargs)

        log.debug(self.args.get('gene_file'))
        cherrypy.log('gene_file:%s'%str(self.args.get('gene_file')))
        if (self.args.get('gene_file') is not None and
            self.args.get('gene_file').file is not None) :
            gene_list = self.args.get('gene_file').file.read()
        else :
            gene_list = self.args.get('gene_list')
        hid = self.args.get('hid')

        top = ''
        if gene_list is None and hid is None :
            top = self.adipo_form()
        else :
            top = self.adipo_results(gene_list,hid)

        if self.args.get('nolog') :
            bottom = ''
        else :
            bottom = self.log(lines=self.args.get('log'))

        template = self.template_env.get_template('adipo_site.html')
        return template.render(top_content=top,bottom_content=bottom)

    def adipo_form(self) :
        log.info('new search')
        top_template = self.template_env.get_template('adipo_site_form.html')
        top_content = top_template.render()
        return top_content

    @check_db_session
    def adipo_results(self,gene_list_str=None,hid=None) :

        cherrypy.log('checking for running thread')
        thread = self.threads[session.id]
        if thread.is_done() :
            top_template = self.template_env.get_template('adipo_site_results.html')
            top_content = top_template.render(**self.templ_d)
        elif thread.is_running() :
            top_content = 'currently running'
            cherrypy.log('thread is running')
        else :
            cherrypy.log('starting thread')
            top_content = 'starting run'
            self.adipo_run(gene_list_str,hid)

        return top_content

    @check_db_session
    def adipo_run(self,gene_list_str=None,hid=None) :

        # gene_list_str comes from the original form
        # hid is a hash of the gene list that is used for recalling
        # gene sets that have already been analyzed
        self.session_d = session_d = {}
        if gene_list_str is not None :

            gene_list_strip = re.sub(r'\s+',' ',gene_list_str)
            log.info('user submitted genes: %s'%gene_list_strip)

            gene_list = gene_list_strip.split()

            # make sure the genes are unique
            gene_list = list(set(gene_list))
            orig_gene_list = gene_list

            # make all genes lower case as they are in the db
            gene_list = [g.lower() for g in gene_list]

            # sort both the original and lower()ed gene lists in parallel,
            # put them back
            both_gene_lists = zip(gene_list,orig_gene_list)
            both_gene_lists.sort()
            gene_list, orig_gene_list = zip(*both_gene_lists)
            self.gene_name_map = dict(zip(gene_list,orig_gene_list)+zip(orig_gene_list,gene_list))

            self.hid = str(hash(''.join(gene_list)))
            session_d['gene_list'] = gene_list
            session_d['orig_gene_list'] = orig_gene_list
            session_d['gene_name_map'] = self.gene_name_map

        elif hid is not None :
            self.hid = hid
            log.debug('session hid passed: %s'%hid)
            self.session_d = session_d = cherrypy.session.get(hid)
            log.debug('session found: %s'%session_d)
            if session_d is None :
                log.error('hid passed but no session found, aborting')
                return 'Error occurred, please try submitting again'

            gene_list = session_d.get('gene_list')
            orig_gene_list = session_d.get('orig_gene_list')
            self.gene_name_map = session_d.get('gene_name_map')
            if gene_list is None :
                log.error('hid is valid but no gene list was found, aborting')
                return 'Error occurred, please try submitting again'

        self.templ_d = templ_d = {}

        templ_d['hid'] = self.hid

        # check for sortby argument from user
        templ_d['sortby'] = self.args.get('sortby')

        # genes found and missing
        templ_d['found'] = []
        templ_d['missing'] = []

        self.region_sets = []
        for gene, orig_gene in zip(gene_list, orig_gene_list) :
            region_set = self.db_session.query(db.RegionSet).filter(db.RegionSet.name==gene).first()
            if region_set is None :
                log.error('could not find record for gene %s'%gene)
                templ_d['missing'].append(orig_gene)
                continue
            else :
                if len(region_set.regions) == 0 :
                    log.error('found record for gene %s, but no region, this is wrong'%gene)
                    templ_d['missing'].append(orig_gene)
                    continue
                self.region_sets.append(region_set)
                templ_d['found'].append(orig_gene)

        # some things might be cached from the last submission
        # get them if they exist
        """
        if 'diff_expr' in session_d :
            diff_expr_d = session_d.get('diff_expr')
        else :
            diff_expr_d = self.diff_expr(gene_list)
        templ_d.update(diff_expr_d)
        """
        templ_d['diff_expr_cnts'] = {}

        if 'motif_enrich' in session_d :
            motif_d = session_d.get('motif_enrich')
        else :
            self.thread.run_thread(self.motif_enrich,gene_list)
            #motif_d = self.motif_enrich(gene_list)

        # store everything into the session
        cherrypy.session[self.hid] = session_d
        #log.debug('stored hid %s in session %s'%(self.hid,cherrypy.session[self.hid]))


    @check_db_session
    def diff_expr(self,gene_list) :

        d = {}

        #conditions = [r.name for r in self.db_session.query(db.Condition.name).all()]
        conditions = ('7-2_bg','7-2_nobg')
        header = ['Gene']+conditions
        cond_diff_exp_cnts = defaultdict(int)
        cond_diff_exp = defaultdict(dict)

        for region_set in self.region_sets :

            orig_gene_name = self.gene_name_map[region_set.name]

            for region in region_set.regions :
                log.debug('Diff Expr db.Region: %s'%region.name)
                for data_rec in region.region_data :
                    log.debug('Diff Expr db.Region Data: %s (%s)'%(data_rec.value,data_rec.data_type.name))
                    if data_rec.data_type.name == 'log2 expression fold change' :
                        if float(data_rec.meta2) < self.args.get('diff_exp_pval') and \
                           (data_rec.value < self.args.get('diff_exp_low') or 
                            data_rec.value > self.args.get('diff_exp_hi')) : # meta2 is qvalue
                            cond_diff_exp[region.name]['name'] = orig_gene_name
                            cond_diff_exp[region.name][data_rec.condition.name] = data_rec.value
                            cond_diff_exp_cnts[data_rec.condition.name] += 1

#            # condition data
#            for condition in conditions :
#                cond_data = [reg_data for reg_data in region.region_data if reg_data.condition.name == condition]
#                for r in cond_data :
#                    if r.data_type.name == 'log2 expression fold change' : #                        data[-1].append(r.value)

        #pass_genes = [r[0] for r in data]
        #data_mat = np.array([r[1:] for r in data])
        #heatmap_fn = 'images/heatmap.png'
        #templ_d['heatmap_fn'] = heatmap_fn

        # make the log2 fold change heatmap
        #self.make_heatmap(data_mat,conditions,pass_genes,heatmap_fn,'log2 fold change')

        d['diff_expr_cnts'] = cond_diff_exp_cnts
        d['diff_expr'] = cond_diff_exp
        self.session_d['diff_expr'] = d

        return d

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

    @check_db_session
    def motif_enrich(self,gene_list) :

        d = {}

        # get motif scores
        cherrypy.log('getting motif scores')

        # sqlite has a limitation of 999 SQL variables, need to split up
        # this query
        batch_size = 500
        gene_list_batches = [gene_list[i:i+500] for i in xrange(0,len(gene_list),500)]
        cherrypy.log('gene_list_batches len()s: %s'%str([len(b) for b in gene_list_batches]))
        scores = []
        for batch in gene_list_batches :
            score_q = (self.db_session.query(db.Region,db.SeqData)
                               .join((db.RegionSet,db.Region.region_sets))
                               .join(db.SeqData)
                               .filter(db.SeqData.seq_type.has(db.SeqType.name=='motif scores'))
                               .filter(db.RegionSet.name.in_(batch))
                     )
            scores.extend(score_q.all())

        cherrypy.log('scores: %s %d'%(str(scores[:100]),len(scores)))
        score_mat = []
        gene_names = []
        condition_scores = defaultdict(list)
        cherrypy.log('loading motif scores')
        for region, seqdata in scores :
            condition = seqdata.condition.name
            motif_scores = cPickle.loads(seqdata.value)
            condition_scores[condition].append(motif_scores)
            score_mat.append(motif_scores)
        cherrypy.log('found conditions: %s'%condition_scores.keys())
        score_mat = np.array(score_mat)

        # compare motif scores of requested genes to all hypersensitive regions
        # in the dataset

        cherrypy.log('loading motif background')
        motif_bg_fn = resource_filename('adipo_sight','data/background_motif_scores.npy')

        all_scores = np.load(motif_bg_fn)

        # pick n random background sequences, but seed so the same indices are picked for every
        # different # of input DHS regions
        random.seed("jo mama")
        bg_inds = random.sample(xrange(all_scores.shape[0]),min(score_mat.shape[0]*3,all_scores.shape[0]))
        all_scores = all_scores[bg_inds,:]
        cherrypy.log('score_mat.shape = %s'%str(score_mat.shape))
        cherrypy.log('all_scores.shape = %s'%str(all_scores.shape))
        cherrypy.log('done loading motif background')

        """
        if cherrypy.session.get('all_scores_mat') is not None :
            cherrypy.log('loading from session cache')
            all_scores = cherrypy.session.get('all_scores_mat')
        else :

            # would be much faster to cache this
            cherrypy.log('running query')
            #WARNING: query only selects 7-2_nobg!! hopefully this will be ok
            all_genes_q = (self.db_session.query(db.Region,db.SeqData)
                         .join((db.RegionSet,db.Region.region_sets))
                         .join(db.SeqData)
                         .filter(db.SeqData.seq_type.has(db.SeqType.name=='motif scores'))
                         #.filter(db.SeqData.condition_id == 7)
                        )
            cherrypy.log('done running query')

            #cherrypy.log('# all genes: %d'%all_genes_q.count())
            all_genes = []
            num_genes = all_genes_q.count()
            cherrypy.log('picking genes + loading motif scores')
            i = 0
            all_scores = np.zeros((num_genes,score_mat.shape[1]))
            for region, seqdata in all_genes_q.yield_per(1000) :
                if i % 1000 == 0 :
                    cherrypy.log('loaded %d'%i)
                if i >= num_genes :
                    break
                scores = cPickle.loads(seqdata.value)
                all_scores[i] = scores
                i += 1

            cherrypy.log('loaded %d gene sequences'%len(all_genes))

            #for i, (region, seqdata) in enumerate(all_genes) :
            #    scores = cPickle.loads(seqdata.value)
            #    all_scores[i] = scores
            cherrypy.session['all_scores_mat'] = all_scores
            np.save('/tmp/all_scores_mat.npy',all_scores)
        """

        # walk through the conditions and compute scores
        sig_scores = defaultdict(dict)
        hs_regions = {}
        enriched_motifs = {}
        for c, scores in condition_scores.items() :

            hs_regions[c] = len(scores)

            # calculate MWW
            score_mat = np.array(scores)
            pvals = mww_multiprocess(score_mat.T,all_scores.T,True)
            log.debug('motif scores for condition: %s'%c)
            log.debug('pvals: %s'%str(pvals.shape))

            motif_name_fn = resource_filename('adipo_sight','data/motif_names.txt')
            motif_names = np.array(open(motif_name_fn).readlines())

            motif_cluster_fn = resource_filename('adipo_sight','data/motif_clusters.txt')
            motif_cluster_map = dict((i,int(m)) for i,m in enumerate(open(motif_cluster_fn)))

            thresh = pvals < self.args.get('diff_hyp_pval')
            thresh_inds = np.where(thresh)[0]
            thresh_names, thresh_pvals = motif_names[thresh], pvals[thresh]
            thresh_imgs = np.array(['images/motif_logos/%03d_motif.png'%i for i in thresh_inds])
            cluster_set = set()
            for i,n,p in zip(thresh_inds,thresh_names,thresh_pvals) :
                cluster_i = motif_cluster_map[i]
                sig_scores[cluster_i].setdefault('name',set()).add(n.strip())
                sig_scores[cluster_i][c] = min(sig_scores[cluster_i].get(c,1.),p)
                cluster_set.add(cluster_i)

            enriched_motifs[c] = len(cluster_set)

        d['motifs'] = dict(sig_scores)
        d['hs_regions'] = hs_regions
        d['enriched_motifs'] = enriched_motifs

        self.session_d['motif_enrich'] = d

        self.templ_d.update(d)

    def make_boxplot(self,vectors,fn,title=None) :
        fig = mp.figure()
        ax = fig.gca()
        ax.boxplot(vectors,sym='')
        ax.set_title(title)
        fig.savefig(fn)

    def log(self,lines=100) :
        try :
            lines = int(lines)
        except TypeError :
            lines = 100
        return log.get_log_content(lines=lines)

    @cherrypy.expose
    @check_db_session
    def ref(self) :

        template = self.template_env.get_template('ref.html')
        return template.render(conditions = self.db_session.query(db.Condition),
                               region_types = self.db_session.query(db.RegionType),
                               seq_types = self.db_session.query(db.SeqType)
                              )
if __name__ == "__main__" :

    import adipo_sight
    import pkg_resources
    current_dir = os.path.dirname(os.path.abspath(__file__))
    package_dir = pkg_resources.resource_filename(adipo_sight.__name__,'data')
    conf = {'/': {'tools.sessions.on': True,
                  'tools.sessions.name': 'adipo_sight',
                  'tools.sessions.storage_type': "file",
                  'tools.sessions.storage_path': "/tmp",
                  'tools.sessions.timeout': 60,
                  },
            '/images': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(package_dir, 'images'),
                        'tools.staticdir.content_types': {'png': 'image/png'}
                       }
           }
    cherrypy.engine.timeout_monitor.unsubscribe()
    cherrypy.config.update({ 'server.socket_host': '127.0.0.1', # localhost
                             #'server.socket_host': '18.68.8.98', # for forwarding from web server
                             'server.socket_port': 8888
                           })
    cherrypy.quickstart(AdipoSite(), '/', config=conf)
