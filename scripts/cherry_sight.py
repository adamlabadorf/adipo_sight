#!/usr/bin/env python

import cPickle
import logging
import os
import pprint
import random
import re
import sys
import threading

from collections import defaultdict
from itertools import cycle

import cherrypy
#import matplotlib.pyplot as mp
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

class MySessions(object) :

    def __init__(self,storage_dir='/tmp',timeout=None) :
        self.storage_dir = storage_dir
        self.session_fns = {}

    def get(self,hid) :
        if hid in self.session_fns and os.path.exists(self.session_fns[hid]) :
            d = cPickle.load(open(self.session_fns[hid]))
        else :
            d = {}
        return d
    
    def save(self,hid,obj) :
        fn = self.session_fns.get(hid,os.path.join(self.storage_dir,'session%s'%hid))
        self.session_fns[hid] = fn
        cPickle.dump(obj,open(fn,'w'))


class AdipoThread(object) :

    def __init__(self) :
        self.thread = None

    def run_thread(self,thread) :

        if self.is_running() :
            raise Exception('Only one thread may be run at a time')

        self.thread = thread
        self.thread.start()

    def is_running(self) :
        return self.thread is not None and self.thread.is_alive()

    def is_done(self) :
        return self.thread is not None and not self.thread.is_alive()

class EnrichThread(threading.Thread) :

    def __init__(self,session,db_session,args) :
        threading.Thread.__init__(self)
        self.session_d = session
        self.db_session = db_session
        self.args = args

    def run(self) :

        session_d = self.session_d

        if session_d is None :
            log.error('hid passed but no session found, aborting')
            return 'Error occurred, please try submitting again'

        gene_list = session_d['gene_list']

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
                               .join(db.RegionMembership)
                               .join(db.RegionSet)
                               .join(db.SeqData)
                               .filter(db.SeqData.seq_type.has(db.SeqType.name=='motif scores'))
                               .filter(db.RegionSet.name.in_(batch))
                               .filter(db.RegionMembership.dist_to_feature.between(
                                        -int(self.args['upstream']),
                                        int(self.args['downstream'])
                                      ))
                     )
            scores.extend(score_q.all())

        score_mat = []
        gene_names = set()
        condition_scores = defaultdict(list)
        cherrypy.log('loading motif scores')
        for region, seqdata in scores :
            for region_membership in region.membership :
                region_set = region_membership.region_set
                gene_names.add(region_set.name)
            condition = seqdata.condition.name
            motif_scores = cPickle.loads(seqdata.value)
            condition_scores[condition].append(motif_scores)
            score_mat.append(motif_scores)
        cherrypy.log('found conditions: %s'%condition_scores.keys())
        score_mat = np.array(score_mat)

        session_d['found'] = gene_names
        session_d['missing'] = [g for g in gene_list if g not in gene_names]

        # compare motif scores of requested genes to all hypersensitive regions
        # in the dataset

        cherrypy.log('loading motif background')
        # pick n random background sequences, but seed so the same indices are picked for every
        # different # of input DHS regions
        random.seed("jo mama")

        # walk through the conditions and compute scores
        sig_scores = defaultdict(dict)
        hs_regions = {}
        enriched_motifs = {}
        for c, scores in condition_scores.items() :

            score_mat = np.array(scores)

            hs_regions[c] = len(scores)

            # get the background out for this condition
            motif_bg_fn = resource_filename('adipo_sight','data/%s_hypersensitive_peaks_bg_motif_scores.npy'%c)

            all_scores = np.load(motif_bg_fn).T
            bg_inds = random.sample(xrange(all_scores.shape[0]),min(score_mat.shape[0],all_scores.shape[0]))
            all_scores = all_scores[bg_inds,:]

            cherrypy.log('score_mat.shape = %s'%str(score_mat.shape))
            cherrypy.log('all_scores.shape = %s'%str(all_scores.shape))
            cherrypy.log('done loading motif background')

            # calculate MWW
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

        session_d['motif_enrich'] = d

class AdipoSite:

    def __init__(self) :

        # get a database connection and look for genes
        self.db_name = db_name = resource_filename('adipo_sight','data/adipo_sight.db')
        cherrypy.log('loading db session from %s'%self.db_name)
        self.db_session = db.get_session(db_name)

        loader = PackageLoader('adipo_sight','data/tmpl')
        self.template_env = Environment(loader=loader)
        self.template_env.filters['dictofdictsort'] = dictofdictsort
        self.template_env.globals['repr'] = repr
        self.template_env.globals['zip'] = zip

        self.args = {}

        self.threads = defaultdict(AdipoThread)

        self.sessions = MySessions()

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
                     'diff_hyp_pval':0.05,
                     'upstream':10000,
                     'downstream':10000}
        self.args.update(kwargs)

        return self.args

    @cherrypy.expose
    def index(self,**kwargs) :

        cherrypy.log('loading index')

        log.info('new search')
        top_template = self.template_env.get_template('adipo_site_form.html')
        top = top_template.render()

        bottom = ''
        #if self.args.get('nolog') :
        #    bottom = ''
        #else :
        #    bottom = self.log(lines=self.args.get('log'))

        template = self.template_env.get_template('adipo_site.html')
        return template.render(top_content=top,bottom_content=bottom)


    @cherrypy.expose
    def submit(self,gene_list,gene_file=None,**kwargs) :

        if gene_list == 'example' :
            example_fn = resource_filename('adipo_sight','data/p65_bound_gene_list.txt')
            gene_list = open(example_fn).read()

        self.parse_kwargs(kwargs)

        log.debug(self.args.get('gene_file'))
        cherrypy.log('gene_file:%s'%str(gene_file))
        if (gene_file is not None and
            gene_file.file is not None) :
            gene_list_str = gene_file.file.read()
        else :
            gene_list_str = gene_list

        gene_list_strip = re.sub(r'\s+',' ',gene_list_str)
        cherrypy.log('user submitted genes: %s'%gene_list_strip)

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

        hid = str(hash(''.join(gene_list)))
        session_d = self.sessions.get(hid)
        session_d['hid'] = hid
        session_d['gene_list'] = gene_list
        session_d['orig_gene_list'] = orig_gene_list
        session_d['gene_name_map'] = self.gene_name_map

        self.gene_list = gene_list
        self.orig_gene_list = orig_gene_list

        self.sessions.save(hid,session_d)

        # now that we've set everything, redirect to the processing page
        raise cherrypy.HTTPRedirect("http://fraenkel.mit.edu/adipo_sight/processing?hid=%s"%hid,status=303)


    @cherrypy.expose
    def processing(self,hid,*args,**kwargs) :

        hid = str(hid)
        cherrypy.log('processing hid %s (%s)'%(hid,type(hid)))
        session_d = self.sessions.get(hid)
        if session_d is None :
            log.error('hid passed but no session found, aborting')
            return 'Error occurred, please try submitting again'

        # we use a thread to process in the bg
        thread = self.threads.get(hid)

        # processing hasn't begun yet
        if thread is None :

            gene_list = session_d.get('gene_list')
            orig_gene_list = session_d.get('orig_gene_list')
            self.gene_name_map = session_d.get('gene_name_map')
            if gene_list is None :
                log.error('hid is valid but no gene list was found, aborting')
                return 'Error occurred, please try submitting again'

            thread = AdipoThread()
            enrich_thread = EnrichThread(session_d,self.db_session,self.args)
            thread.run_thread(enrich_thread)
            self.threads[hid] = thread

        if thread.is_running() :

            template = self.template_env.get_template('adipo_site.html')
            top = ('Request is processing, please wait.<br/>This page refreshes automatically every 5 seconds.<br/>'
                   '<meta http-equiv="refresh" content="5; URL=processing?hid=%s"</meta>')%hid
            return template.render(top_content=top,bottom_content='')

        elif thread.is_done() :

            raise cherrypy.HTTPRedirect("http://fraenkel.mit.edu/adipo_sight/results?hid=%s"%hid,status=303)


    @cherrypy.expose
    def results(self,hid,*args,**kwargs) :

        args = self.parse_kwargs(kwargs)

        session_d = self.sessions.get(hid)
        if session_d is None :
            log.error('hid passed but no session found, aborting')
            return 'Error occurred, please try submitting again'

        thread = self.threads[hid]
        session_d.update(thread.thread.session_d)
        session_d.update(thread.thread.session_d['motif_enrich'])

        templ_d = {}
        templ_d.update(session_d)

        # check for sortby argument from user
        templ_d['sortby'] = str(args.get('sortby'))

        top_template = self.template_env.get_template('adipo_site_results.html')
        top_content = top_template.render(**templ_d)

        return top_content


    def log(self,lines=100) :
        try :
            lines = int(lines)
        except TypeError :
            lines = 100
        return log.get_log_content(lines=lines)

    @cherrypy.expose
    def motif_logos(self,motifs='') :

        motif_name_fn = resource_filename('adipo_sight','data/motif_names.txt')
        motif_names = open(motif_name_fn).readlines()

        motif_logos = [(n,'%03d_motif.png'%i) for i,n in enumerate(motif_names)]
        cherrypy.log(','.join(motif_names[:10]))

        if len(motifs) != 0 :
            cherrypy.log('got motifs: %s'%motifs)
            motifs = [str(m).strip() for m in motifs.split(',')]
            cherrypy.log('split motifs: %s'%motifs)
            motif_logos = [l for l in motif_logos if l[0].strip() in motifs]

        template = self.template_env.get_template('motif_logos.html')
        return template.render(motif_logos=motif_logos,motifs=motifs)

    @cherrypy.expose
    def query(self,stmt='') :
        db_session = db.get_session(self.db_name)
        resp = ''
        fields = None
        if len(stmt) > 0 :
            if any([k in stmt.lower() for k in ('update','insert','delete')]) :
                resp = 'illegal SQL statment'
            else :
                try :
                    fields = re.search('select(.*?)from',stmt.lower()).group(1)
                    fields = [f.strip() for f in fields.split(',')]
                    recs = db_session.query(*fields).from_statement(stmt)
                    resp = [dict((f,getattr(r,f)) for f in fields) for r in recs]
                except Exception, e :
                    cherrypy.log('exception: %s'%e)
                    resp = str(e.args)

        template = self.template_env.get_template('query.html')
        return template.render(stmt=stmt,resp=resp,fields=fields)

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
    print package_dir
    conf = {'/': {'tools.sessions.on': True,
                  'tools.sessions.name': 'adipo_sight',
                  'tools.sessions.storage_type': "file",
                  'tools.sessions.storage_path': "/tmp",
                  'tools.sessions.timeout': 60,
                  },
            '/images': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(package_dir, 'images'),
                        'tools.staticdir.content_types': {'png': 'image/png'}
                       },
            '/css': {'tools.staticdir.on': True,
                     'tools.staticdir.dir': os.path.join(package_dir, 'css'),
                     'tools.staticdir.content_types': {'css': 'text/css'}
                    }
           }
    cherrypy.config.update({ 'server.socket_host': '127.0.0.1', # localhost
                             #'server.socket_host': '18.68.8.98', # for forwarding from web server
                             'server.socket_port': 8888
                           })
    cherrypy.quickstart(AdipoSite(), '/', config=conf)
