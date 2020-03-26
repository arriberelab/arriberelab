"""
Joshua Arribere, March 26, 2020

Script to look at different ways of making a python dict and seeing
    which is most efficient.

Input: nothing

Output: print-to-screen of dict's memory fp

run as python testDictStructures.py
"""
import sys
from logJosh import Tee
from random import randrange
import os, psutil

def method1(chrs,chrSize,N):
    aa={}
    for ii in range(chrs):
        Chr='chr%s'%(ii)
        if Chr not in aa:
            aa[Chr]={}
        for jj in range(N):
            kk=randrange(0,chrSize)
            if kk not in aa[Chr]:
                aa[Chr][kk]=[]
            aa[Chr][kk].append('Somestring:%s:S'%(randrange(10**5)))
    return aa

def method2(chrs,chrSize,N):
    #very similar to method1 in overall usage
    aa={}
    for ii in range(chrs):
        Chr='chr%s'%(ii)
        if ii not in aa:
            aa[ii]={}
        for jj in range(N):
            kk=randrange(0,chrSize)
            if kk not in aa[ii]:
                aa[ii][kk]=[]
            aa[ii][kk].append('Somestring:%s:S'%(randrange(10**5)))
    return aa

def method3(chrs,chrSize,N):
    #More than either method1 or 2
    aa={}
    for ii in range(chrs):
        for jj in range(N):
            kk=randrange(0,chrSize)
            if (ii,kk) not in aa:
                aa[(ii,kk)]=[]
            aa[(ii,kk)].append('Somestring:%s:S'%(randrange(10**5)))
    return aa

def method4(chrs,chrSize,N):
    #More than either method1 or 2, but not as much as 3
    aa={}
    for ii in range(chrs):
        for jj in range(N):
            kk=randrange(0,chrSize)
            key='chr%s_%s'%(ii,kk)
            if key not in aa:
                aa[key]=[]
            aa[key].append('Somestring:%s:S'%(randrange(10**5)))
    return aa

def method5(chrs,chrSize,N):
    #similar in usage to method4
    aa={}
    for ii in range(chrs):
        Chr='chr%s'%(ii)
        if Chr not in aa:
            aa[Chr]={}
        for jj in range(N):
            kk=randrange(0,chrSize)
            key=str(kk)
            if key not in aa[Chr]:
                aa[Chr][key]=[]
            aa[Chr][key].append('Somestring:%s:S'%(randrange(10**5)))
    return aa

def method6(chrs,chrSize,N):
    #a little over half the mem usage of method1
    aa={}
    for ii in range(chrs):
        Chr='chr%s'%(ii)
        if Chr not in aa:
            aa[Chr]={}
        for jj in range(N):
            kk=randrange(0,chrSize)
            if kk not in aa[Chr]:
                aa[Chr][kk]=''
            aa[Chr][kk]+='Somestring:%s:S'%(randrange(10**5))
    return aa

def method7(chrs,chrSize,N):
    #between methods 1 and 6
    aa={}
    for ii in range(chrs):
        Chr='chr%s'%(ii)
        if Chr not in aa:
            aa[intern(Chr)]={}
        for jj in range(N):
            kk=randrange(0,chrSize)
            kk=str(kk)
            if kk not in aa[intern(Chr)]:
                aa[intern(Chr)][intern(kk)]=''
            aa[intern(Chr)][intern(kk)]+='Somestring:%s:S'%(randrange(10**5))
    return aa

def method8(chrs,chrSize,N):
    #essentially identical to method6
    aa={}
    for ii in range(chrs):
        Chr='chr%s'%(ii)
        if Chr not in aa:
            aa[intern(Chr)]={}
        for jj in range(N):
            kk=randrange(0,chrSize)
            if kk not in aa[intern(Chr)]:
                aa[intern(Chr)][kk]=''
            aa[intern(Chr)][kk]+='Somestring:%s:S'%(randrange(10**5))
    return aa


def main(args):
    #first compute as nested dict
    chrs=20
    chrSize=10**9
    N=10**5
    #
    aa=method1(chrs,chrSize,N)
    #
    #bb=method2(chrs,chrSize,N)
    #
    #cc=method3(chrs,chrSize,N)
    #
    #dd=method4(chrs,chrSize,N)
    #
    #ee=method5(chrs,chrSize,N)
    #
    #ff=method6(chrs,chrSize,N)
    #
    #gg=method7(chrs,chrSize,N)
    #
    #hh=method8(chrs,chrSize,N)
    #
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss/(10**6))


if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
