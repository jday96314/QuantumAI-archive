"""Fixed-word rank-C refit control with explicit spectral-work budget."""
import argparse
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--budget',type=int,default=6000);a=p.parse_args();z,u,w,tr,words=load(a.init);R=w.shape[0]
    start=time.time();m=Refit(u);initial=physical_error(u,w);cycles=(a.budget+R-1)//R
    w,err,hist=m.fit(w,cycles)
    save(HERE/f'candidates/{a.tag}.npz',u,w,tr,words,err2=err)
    report(HERE/f'results/{a.tag}.json',dict(initial=initial,final=err,requested_block_eigensolves=a.budget,
           actual_cycles=len(hist)-1,upper_bound_block_eigensolves=cycles*R,history=hist,seconds=time.time()-start))
if __name__=='__main__':main()
