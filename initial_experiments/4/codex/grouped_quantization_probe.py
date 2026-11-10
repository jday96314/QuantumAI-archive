"""Measure direct quantization error for grouped dictionaries before circuit work.

This probe preserves each reference weight array. Later fixed-dictionary PSD
refits can improve the error; none of these measurements is native acceptance.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from angle_gauge_probe import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors


def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--bits',type=int,nargs='+',default=[7,8,9,10,12]);p.add_argument('--stem',required=True);a=p.parse_args();torch.set_num_threads(2);begin=time.time()
 z=dict(np.load(a.input));sp=np.load('../../3/parameters/reference_spectrum.npz');d=dict(device='cuda',dtype=torch.float64)
 G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d);ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];records=[]
 for bits in a.bits:
  words=encode(z['U'].reshape(-1,76),z['tree_edges'],bits,True);U=reconstruct_u_vectors(words/2**(bits-1),z['tree_edges']).reshape(z['U'].shape)
  P=torch.tensor(U[:,:,ii]*U[:,:,jj]*fac,**d);X=(torch.tensor(z['w'],**d).transpose(1,2)@P).reshape(-1,len(ii));er=float((X.T@X-G).square().sum().sqrt())
  path=f'parameters/{a.stem}_q{bits}_seed.npz';out=dict(z);out.update(U=U,angle_words=words,brot=bits,err2=er);np.savez(path,**out)
  row=dict(bits=bits,error_2=er,output=path,R=U.shape[0],B=U.shape[1],C=z['w'].shape[2]);records.append(row);print(json.dumps(row),flush=True)
 Path(f'results/{a.stem}_quantization_probe.json').write_text(json.dumps(dict(scope='unrefitted direct angle quantization; no native acceptance',input=a.input,records=records,seconds=time.time()-begin),indent=2)+'\n')
if __name__=='__main__':main()
