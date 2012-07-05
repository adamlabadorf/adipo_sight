""" This code written by Christopher Ng, adapted by Adam Labadorf """

from __future__ import division
import sys, math, cPickle,os
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from optparse import OptionParser
from collections import defaultdict

def calc_bh(p_values, num_total_tests):
    """
    Calculates the Benjamini-Hochberg correction for multiple hypothesis
    testing from a list of p-values *sorted in ascending order*.

    """
    BH=[]
    prev_bh_value = 0
    for i, p_value in enumerate(p_values):
        bh_value = p_value * num_total_tests / (i + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1)

        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        #yield bh_value
        BH.append(bh_value)
    return BH 


def mww(A,B,BH=False) :
    """Perform Mann-Whitney U test on matrices A and B.  A pvalue is produced
    for every row of the matrices (A and B must have the same number of
    rows.)  BH == True does Benjamini-Hochberg multiple hypothesis correction
    under independence assumption.  Returns p-values in the order of the
    original matrices."""

    n_motifs, n1 = A.shape
    n_motifs, n2 = B.shape

    Z = np.zeros((n_motifs,3))
    for i in np.arange(n_motifs):

        m_a, m_b = A[i,:], B[i,:]
        #n1,n2=len(A),len(B)
        M_a_label=[(x,True) for x in m_a]
        M_b_label=[(x,False) for x in m_b]
        M=M_a_label+M_b_label
        M.sort()
        M_rev=M[::-1]
        R1,R2=0,0
        
        M_dict=defaultdict(list)
        r=1
        last_score=0
        inc=1
        for score,label in M_rev:
            M_dict[r].append(label)
            if last_score!=score: 
                r+=inc
                inc=1
            else: inc+=1
            last_score=score
        tied_ranks=len(M_dict.keys())
        TIES=[]
        for k,v in M_dict.items():
            ties=len(v)
            TIES.append(ties)
            mean_rank=k+(ties-1)/2
            for val in v:
                if val: R1+=mean_rank
                else: R2+=mean_rank

        U1=R1-n1*(n1+1)/2
        U2=R2-n2*(n2+2)/2
        AUC1=U1/(n1*n2)
        AUC2=U2/(n1*n2)
        mu_U=n1*n2/2

        t_k=0
        for t in TIES: t_k+=(t**3-t)

        n=n1+n2
        sigma_U=math.sqrt(n1*n2/n/(n-1))*math.sqrt((n**3-n-t_k)/12)
        z1=(U1-mu_U)/sigma_U
        z2=(U2-mu_U)/sigma_U
        pval=scipy.stats.norm.cdf(z1)

        Z[i,:] = i,z1,pval #,'%i'%i,id,'%f'%z1,'%f'%z2,'%i'%U1,'%i'%U2,'%f'%AUC1,'%f'%AUC2])

    #BH correction
    if BH :

        # sort in ascending pval order
        Z = Z[Z[:,2].argsort(),:]

        Z[:,2] = calc_bh(Z[:,2],n_motifs)

        # sort back
        Z = Z[Z[:,0].argsort(),:]

    # only return the p-values
    Z = Z[:,2]
    return Z


def main():
    
    from numpy.random import normal

    n, m = 3,20
    a = np.zeros((n,m))
    a[0,:] = normal(5,1,(m))
    a[1,:] = normal(3,1,(m))
    a[2,:] = normal(10,1,(m))

    m = 100
    b = np.zeros((n,m))
    b[0,:] = normal(3,1,(m))
    b[1,:] = normal(3,1,(m))
    b[2,:] = normal(3,1,(m))

    print 'a means: ',a.mean(1)
    print 'b means: ',b.mean(1)
    
    pvalues = mww(a,b,False)
    print 'non-BH corrected:',pvalues

    pvalues = mww(a,b,True)
    print 'BH corrected:',pvalues


if __name__=='__main__':

    main()
