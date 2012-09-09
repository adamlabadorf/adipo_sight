from __future__ import division
import sys,math
from TAMO import MotifTools
pypath=['/nfs/data/cwng/python_code/cvEM.64/orig/cvEM.64']
for x in range(0,len(pypath)):  sys.path.append(pypath[x])
#from TAMO import MotifTools
import numpy as np
from optparse import OptionParser
#import MDsupport
import EM
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from TAMO.seq import Fasta

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--genome", dest="genome", default='mm9')
    parser.add_option("--outfile", dest="outfile")
    parser.add_option("--motif", dest="motif",default="/nfs/vendata/cwng/TRANSFAC_2011.3/TAMO_filt/vertebrates_all_ic8.tamo")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    fsa=args[0]
    motif_matrix(fsa,options.motif,options.outfile,genome=options.genome)

def motif_matrix(fsa,motif,outfile,genome='mm9'):
    if genome=='hg18': markov="/nfs/genomes/human_gp_mar_06/hg18_promoters_3000_1000.markov"
    else: markov="/nfs/data/cwng/chipseq/hypotheses/Mouse.markov"

    #Load motif and background adjust PSSM
    m=MotifTools.load(motif)
    EM.loadMarkovBackground(markov)
    bg = EM.theMarkovBackground.zeroth()
    F=Fasta.load(fsa,key_func=lambda x:x)
    seqs=F.values()
    n_seqs=len(seqs)
    n_motifs=len(m)
    SCORES=np.zeros((n_motifs,n_seqs),dtype='float')
    #SHIFTS=np.zeros((n_motifs,n_seqs))
    
    #out=open(outfile,'w')
    for i,M in enumerate(m):
        ll = M.logP
        EM.loadMarkovBackground(markov)
        bg = EM.theMarkovBackground.zeroth()
        for pos in ll:
            for letter in pos.keys():
                pos[letter] = pos[letter] - math.log(bg[letter])/math.log(2.0)
        AM = MotifTools.Motif_from_ll(ll)
        #adj_model = MotifTools.Motif_from_ll(ll)
        #adj_model.source = M.source
        #pssm = MDsupport.Motif2c_PSSM(adj_model)
        #w=pssm.width

        #shift=[]
        #scores=[]
        mi,ma=AM.minscore,AM.maxscore

        #F_m={}
        #Search every seq for given motif above threshold t and print motif centered results
        for j,seq in enumerate(seqs):
            seq_fwd = seq.upper()
            #seq_rev = str(MotifTools.revcomplement(seq_fwd))[::-1]
            #scores_fwd = pssm.score_probe(seq_fwd)
            #scores_rev = pssm.score_probe(seq_rev)
            #max_score=mi
            #max_ind=0
            #for ind,s in enumerate(scores_fwd):
            #    if s> max_score:    
            #        max_score=s
            #        max_ind=ind
            #        strand='+'
            #for ind,s in enumerate(scores_rev):
            #    if s> max_score:
            #        max_score=s
            #        max_ind=ind
            #        strand='-'
            max_score=AM.bestscore(seq_fwd)
            mscore=(max_score-mi)/(ma-mi)
            #orig=len(seq_fwd)/2
            #bind=max_ind+w//2
            #d=abs(orig-bind)
            SCORES[i,j]=mscore
            #SHIFTS[i,j]=d
            #out.write('%1.3f\t'%mscore)
        #out.write('\n')
    #out.close()
    #del F
    np.savetxt(outfile,SCORES,fmt='%1.3f')

if __name__=='__main__': main()
