"""Prototype in-place Bezout(E,O) multiplication using the upstream Euclid tape.

Adapted from Apache-2.0 andigu/uncompute, pinned revision documented in README.
Source polynomials are restored; one payload vector is transformed in place.
The current prototype assumes ordinary degree-one Euclidean quotients.
"""
import inspect,json,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
ROOT=Path(__file__).resolve().parent

def build(full=False):
 m=load();oldswap=m.Em.swap_fn
 m.Em.swap_fn=lambda self,control,pairs:oldswap(self,control,tuple(pairs))
 b=m.Builder();EW=[b.qubits(12) for _ in range(18)];OW=[b.qubits(12) for _ in range(18)];X=[b.qubits(12) for _ in range(17)]
 V=[b.qubits(12) for _ in range(18)] if full else []
 tag='sylvester' if full else 'bezout'
 A=m.Alloc(b);sA,sB,sT=b.qubits(3);E=m.Em(b,sA,sB,sT,A)
 eps=A.open('eps',6);pre=A.open('pre',1)[0];JF=m.Journal();JPRE=m.Journal();records=[];counts=[]
 def scale(word,coefficient):
  with m.inv_word(E,coefficient) as iv,m.opened(E,'bproduct',12) as prod:m.mon_slot(E,word,coefficient,iv,prod)
 def middle(coeff,low,high):
  n=len(low)+len(high)-1;pts,I=m.vinv(n)
  for a,point in enumerate(pts):
   ca=[m.cmrows(I[j][a]) for j in range(n)]
   cb=[m.cmrows(m.gpow(point,j)) for j in range(len(low))]
   sa=[frozenset(coeff[j+1][k] for j in range(n) for k in m.bitlist(ca[j][bit],12)) for bit in range(12)]
   sb=[frozenset(low[-1-j][k] for j in range(len(low)) for k in m.bitlist(cb[j][bit],12)) for bit in range(12)]
   powers=[m.gpow(point,i) for i in range(len(high))]
   def outputs(tb):
    return [high[i][k] for i,power in enumerate(powers) for k in m.bitlist(m.gmul(power,m.tmask(tb)),12)]
   terms=[];m.bil_terms(terms,m.KT,sa,sb,tws_of=outputs);E.emit_mac(terms)
 def triangular(words,inverse=False):
  if len(words)<=1:return
  cut=len(words)//2;lo,hi=words[:cut],words[cut:]
  if inverse:triangular(lo,True);middle(EW,lo,hi);triangular(hi,True)
  else:triangular(hi);middle(EW,lo,hi);triangular(lo)
 if full:
  triangular(V)
  terms=m.polyprod_rt(OW[:17],X,[[V[k+1]] if k<17 else [] for k in range(33)])
  E.emit_mac(terms)
  triangular(X,True)
  E.swaps(None,m.word_pairs(X[:8],list(reversed(X))[:8]))
 def forward_callback(em,R,qh,j,marks):
  d=len(R)-1;records.append((list(R),dict(marks)))
  start=b._stream_totals()[0]
  with E.using(J=m.NOJ):
   scale(X[d],R[-1])
   for i in range(d):E.emit_mac(m.prod_rt(X[d],R[i],X[i]))
   scale(X[d],qh)
  counts.append(b._stream_totals()[0]-start)
 namespace=dict(m.__dict__);namespace['payload_callback']=forward_callback
 def quotient_inplace(em,word,inverse,pivot):
  with m.opened(em,'qlproduct',12) as prod:
   em.J.step(lambda:m.mon_slot(em,word,inverse,pivot,prod),
             lambda:m.mon_slot(em,word,pivot,inverse,prod))
 namespace['quotient_inplace']=quotient_inplace
 source=inspect.getsource(m.psi_forward)
 source=source.replace('def psi_forward(', 'def payload_forward(')
 start=source.index('            with opened(E, "qh", B) as qh, opened(E, "ql", B) as ql:')
 end=source.index('            g = g[16:] + g[:16]',start)
 source=source[:start]+'''            with opened(E, "qh", B) as qh:
                E.mac(prod_rt(qh, g[17], fi))
                if j <= 16:E.xor(prod_rt(g[16], qh, f[16]))
                ql = g[16]
                quotient_inplace(E, ql, fi, f[17])
                mk["q"] = E.J.mark()
                if j <= 16:
                    targets=[[g[j-1+n]] if j-1+n < 16 else [] for n in range(18-j)]
                    E.xor(polyprod_rt([ql, qh], [f[k] for k in range(j,17)], targets))
                E.free(prod_rt(g[17], qh, f[17]))
                E.swaps(None, word_pairs([qh],[g[17]]))
''' +source[end:]
 source=source.replace('recs[j] = {"marks": mk, "flag": flag, "c": c, "q": [ql, qh]}','recs[j] = {"marks": mk, "flag": flag, "c": c, "q": [ql, qh]}\n        payload_callback(E, f[j:], g[1], j, mk)')
 exec(source,namespace)
 with E.using(J=JPRE):m.psi_pre(E,EW,OW,pre)
 with E.using(J=JF):namespace['payload_forward'](E,EW,OW,eps[0],eps[1:])
 # Callback sites are at clean scratch boundaries, so scratch IDs can be reused
 # between source-journal blocks without holding the full historical pool.
 for R,marks in reversed(records):
  d=len(R)-1
  with E.using(J=m.NOJ):
   for i in range(d):E.emit_mac(m.prod_rt(X[i],R[i],X[d]))
   scale(X[d],R[-1])
  JF.reverse(marks['start'],marks['end'])
 first=records[0][1]['start'];JF.reverse(0,first);JPRE.reverse()
 if full:
  triangular(V,True)
  terms=m.polyprod_rt(EW[1:],V[1:],[[] if k<16 else [X[k-16]] for k in range(33)])
  E.emit_mac(terms);triangular(V)
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
 tape=ROOT/f'parameters/{tag}_tape.npy';np.save(tape,rows)
 layout=dict(E=EW,O=OW,X=X,V=V,maxq=b.next_qubit,outcomes=b.next_out,ops=len(rows),T=b._stream_totals()[0],payload_forward_T=sum(counts))
 (ROOT/f'parameters/{tag}_hook.json').write_text(json.dumps(layout,indent=2))
 prefile=ROOT/f'parameters/{tag}.pre.kmx';kb=km.CircuitBuilder()
 for name,words in [('E',EW),('O',OW),('X',X),('V',V)]:
  for i in range(len(words)):kb.create_quantum_register(12,f'{name}{i}')
 kb.finish_circuit().to_file(str(prefile),format='kmx')
 empty=ROOT/f'parameters/{tag}.empty';empty.write_text('');path=ROOT/f'candidates/{tag}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(prefile),str(empty),str(empty),str(path),str(b.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb')
 print(json.dumps({k:v for k,v in layout.items() if k not in ['E','O','X','V']}),flush=True)

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--sylvester',action='store_true');args=parser.parse_args();build(args.sylvester)
