"""Approximate completion chart with explicit X-outcome bookkeeping.

Coincident nonzero active heads and accidental identity affine factors are
unsupported rare charts. Complete native verification measures their cost.
"""
from contextlib import contextmanager
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from extension_division import create as divide_create,phase as divide_phase
from schur_dynamic_fold import dapply,powers
from schur_dynamic_division_contract import phase as dynamic_phase
from schur_dynamic_division_sequential import create as dynamic_create
from schur_fold import axpy,xword
from schur_control import caxpy
from low_workspace_inverse import bind
from fast_dqi import FastBuilder
from schur_affine_base_quotient import create as base_create


@contextmanager
def evaluate(b,E,N,B,kappa,controls):
 m=setup()[0];deg=len(B);assert len(N)==deg+1 and deg>=1;extra=None;nv=N[:2]
 if deg==1:extra=b.alloc_qubits(12);bv=[B[0],extra]
 else:bv=B[:2]
 pp=powers(deg,kappa,m)
 def linear():
  for j in range(2,deg+1):axpy(b,N[j],nv[0],pp[j][0],m);axpy(b,N[j],nv[1],pp[j][1],m)
  for j in range(2,deg):axpy(b,B[j],bv[0],pp[j][0],m);axpy(b,B[j],bv[1],pp[j][1],m)
  xword(b,bv[0],pp[deg][0]);xword(b,bv[1],pp[deg][1])
 linear()
 def correction(flag,delta):
  alt=powers(deg,kappa^delta,m);dif=[(a^c,b^e) for (a,b),(c,e) in zip(pp,alt)]
  for bank,dest,monic in [(N,nv,False),(B,bv,True)]:
   higher=list(range(2,len(bank)))
   if len(higher)==1:
    j=higher[0];rt=[]
    for bit,q in enumerate(bank[j]):
     lo=m.gmul(dif[j][0],1<<bit);hi=m.gmul(dif[j][1],1<<bit);targets=[dest[0][i].id for i in range(12) if lo>>i&1]+[dest[1][i].id for i in range(12) if hi>>i&1]
     if targets:rt.append(('P',(flag,),0,(q.id,),0,tuple(targets)))
    E.mac(rt)
   elif higher:
    tmp=b.alloc_qubits(24)
    def summation():
     for j in higher:axpy(b,bank[j],tmp[:12],dif[j][0],m);axpy(b,bank[j],tmp[12:],dif[j][1],m)
    summation()
    for q,t in zip(tmp,list(dest[0])+list(dest[1])):b.ccx(km.q(flag),q,t)
    summation();b.free(tmp)
   if monic:
    for bit in range(12):
     if dif[deg][0]>>bit&1:b.cx(km.q(flag),dest[0][bit])
     if dif[deg][1]>>bit&1:b.cx(km.q(flag),dest[1][bit])
 for flag,delta in controls:correction(flag,delta)
 yield [q.id for w in nv for q in w],[q.id for w in bv for q in w]
 for flag,delta in reversed(controls):correction(flag,delta)
 linear()
 if extra is not None:b.free(extra)


def affine(b,E,u,t,N,B,inverse=False):
 m=setup()[0];m.inv_word=bind(m);U=[q.id for q in u];T=[q.id for q in t]
 with m.inv_word(E,T) as ti,m.opened(E,'affine_chain_product',12) as temp:
  def scale(back=False):
   for w in N:m.mon_slot(E,[q.id for q in w],ti if back else T,T if back else ti,temp)
  def shear():
   for n,v in zip(N[:-1],B):E.mac(m.prod_rt([q.id for q in n],U,[q.id for q in v]))
   b.cx(u,N[-1])
  if inverse:shear();scale(True)
  else:scale();shear()


def build(d,kappa,deltas,inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());headers=[];rid=0
 def reg(n,name):
  nonlocal rid
  q=b.alloc_qubits(n)
  for x in q:headers.append(f'APPEND_TO_REGISTER q{x.id} r{rid}')
  headers.append(f'REGISTER r{rid} "{name}"');rid+=1;return q
 u=reg(12,'u');t=reg(12,'t');N=[reg(12,f'N{i}') for i in range(d+1)];B=[reg(12,f'B{i}') for i in range(d)];nl=reg(12,'Nlead');bt=reg(12,'Btop');completion=b.alloc_qubits(1)[0];flags=reg(len(deltas),'flags');controls=[(q.id,v) for q,v in zip(flags,deltas)];head_ghost=b.alloc_bits(25);headers.extend(f'APPEND_TO_REGISTER b{v.id} r{rid}' for v in head_ghost);headers.append(f'REGISTER r{rid} \"head_ghost\"');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));NN=N+[nl];BB=B+[bt]
 def correction():
  for bit in range(12):
   if kappa>>bit&1:b.cx(completion,N[-1][bit])
  b.cx(completion,bt[0]);caxpy(m,E,completion.id,[q.id for q in nl],[q.id for q in bt],1)
 if inverse:
  if not d:base_create(b,E,N[0],nl,bt,u,t,kappa,controls)
  else:
   with evaluate(b,E,NN,BB,kappa,controls) as (nv,bv):dynamic_create(b,E,nv,bv,[q.id for w in [u,t] for q in w],kappa,controls)
  b.flip_if_equal(t,0,target=completion);b.cx(completion,t[0]);caxpy(m,E,completion.id,[q.id for q in nl],[q.id for q in u],1);affine(b,E,u,t,NN,BB,True);correction();dapply(b,E,N,BB,kappa,controls,True);b.x(bt[0]);b.cx(N[-1],nl)
 else:
  b.bit_store(head_ghost,False)
  b.flip_if_equal(km.array(list(u)+list(t)),1<<12,target=completion);b.cx(N[-1],nl);b.x(bt[0]);dapply(b,E,N,BB,kappa,controls);correction();affine(b,E,u,t,NN,BB)
  b.hmr(km.array(list(u)+list(t)),head_ghost[:24]);b.free(u);b.free(t)
  b.hmr(completion,head_ghost[24]);b.free(km.array([completion]))
 if inverse:
  for q,bit in zip(list(u)+list(t),head_ghost[:24]):b.push_condition(bit);b.z(q);b.pop_condition()
  b.push_condition(head_ghost[24]);b.z(completion);b.pop_condition();b.bit_store(head_ghost,False)
  b.flip_if_equal(km.array(list(u)+list(t)),1<<12,target=completion)
  reset=b.alloc_bits(1)[0];b.hmr(completion,reset);b.bit_store(reset,False);b.free(km.array([reset]));b.free(km.array([completion]))
 ab.release_recycle('affine_chain_emitter',sc);c=b.finish_circuit();return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(25,c.num_bits)))


