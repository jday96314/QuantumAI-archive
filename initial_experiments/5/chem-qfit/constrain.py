# Adapted from round-3 constrained_factor.py for smaller integer dictionaries; native checks remain separate.
"""Optimize a fixed quantized dictionary's factors inside physical error gates."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=6000);p.add_argument('--lr',type=float,default=1e-5)
p.add_argument('--cap',type=float,default=.019995);p.add_argument('--whitened',action='store_true');p.add_argument('--wb-offset',action='store_true')
p.add_argument('--wb-scale',type=float);p.add_argument('--offset-lr-factor',type=float,default=3.)
p.add_argument('--mix-penalty',type=float,default=0.,help='coefficient of sign mixing in physical squared-error units; divided by cap**2 in the normalized objective')
p.add_argument('--mix-normalization',choices=['physical_error','objective'],default='physical_error',help='physical_error scales the coefficient by cap**-2; objective adds the literal coefficient to the existing loss')
p.add_argument('--optimizer',choices=['adam','lbfgs'],default='adam');a=p.parse_args();torch.set_num_threads(3)
offset_scale=a.wb_scale if a.wb_scale is not None else float(a.wb_offset)
mix_coefficient=a.mix_penalty/a.cap**2 if a.mix_normalization=='physical_error' else a.mix_penalty
HERE=Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3');z=dict(np.load(a.fit));u=z['U'][0];w0=z['w'][0];W0=z['W'][0]
sp=np.load(HERE/'parameters/reference_spectrum.npz');ref=np.load(HERE/'parameters/chem_reference.npz')
ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];P=u[:,ii]*u[:,jj]*fac
ev,V=np.linalg.eigh(P@P.T);root=(V*np.sqrt(ev))@V.T;iroot=(V/np.sqrt(ev))@V.T;Q=iroot@P
G=(sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T;H=Q@G@Q.T
outside=float(np.square(sp['eigenvalues']).sum()-np.square(H).sum())
d=dict(device='cuda',dtype=torch.float64)
rt=torch.tensor(root,**d);irt=torch.tensor(iroot,**d);ht=torch.tensor(H,**d);pt=torch.tensor(P,**d);ft=torch.tensor(fac,**d)
ti=torch.tensor(ii,device='cuda');tj=torch.tensor(jj,device='cuda')
w=torch.tensor(root@w0 if a.whitened else w0,requires_grad=True,**d);off=torch.tensor(W0+offset_scale*w0.sum(0),requires_grad=True,**d)
h0=torch.tensor(ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(76),**d)
constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
def vals():
 ww=irt@w if a.whitened else w
 oo=off-offset_scale*ww.sum(0)
 aa=w if a.whitened else rt@w;e2=outside+(aa@aa.T-ht).square().sum()
 packed=(ww@oo)@pt/ft;hc=torch.zeros((76,76),**d);hc[ti,tj]=packed
 hc=hc+hc.T-torch.diag(hc.diag());h=h0-hc;e=torch.linalg.eigvalsh(h);r=torch.sqrt(e*e+1e-14)
 gap=constant+r.sum()-e.sum()+.5*oo.square().sum()
 eps=1e-7;lam=r.sum()+.25*(torch.sqrt(ww*ww+eps*eps).sum(0)+torch.sqrt((oo+ww.sum(0)).square()+eps*eps)).square().sum()
 score=gap*(2*lam-gap)
 return e2,gap,lam,score
opt=(torch.optim.Adam([{'params':[w],'lr':a.lr},{'params':[off],'lr':a.lr*a.offset_lr_factor}]) if a.optimizer=='adam' else torch.optim.LBFGS([w,off],lr=a.lr,max_iter=12,history_size=25,line_search_fn='strong_wolfe',tolerance_grad=1e-10,tolerance_change=1e-14))
rho=2.;mue=float(z.get('lagrange_error',0.));mug=float(z.get('lagrange_gap',0.));best=float('inf');bestmix=float('inf');start=time.time()
def mix_metrics():
 ww=irt@w if a.whitened else w
 l1sq=ww.abs().sum(0).square().sum();trsq=ww.sum(0).square().sum()
 return l1sq-trsq,l1sq,trsq
def save_fit(path,er):
 actual=irt@w if a.whitened else w;actualoff=off-offset_scale*actual.sum(0)
 out=dict(z)
 for key in ['onebody_angle_words','onebody_brot','error_1_optimized']:out.pop(key,None)
 out.update(w=actual.detach().cpu().numpy()[None,:,:],W=actualoff.detach().cpu().numpy()[None,:],err2=er,lagrange_error=mue,lagrange_gap=mug,mix_penalty=a.mix_penalty,mix_normalization=a.mix_normalization)
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('wb') as f:np.savez(f,**out)
 tmp.replace(path)
with torch.no_grad():
 ee,gg,ll,ss=vals();er=float(ee.clamp(min=0).sqrt());gv=float(gg);sv=float(ss.clamp(min=0).sqrt())
 if er<=.019999 and gv>=3.9151:
  best=sv;bestmix=float(mix_metrics()[0]);save_fit(a.out,er)
  print(json.dumps(dict(step=-1,error_2=er,gap=gv,Lambda=float(ll),lambda_eff=sv,best=best,mixing=bestmix,l1_squared=float(mix_metrics()[1]),mix_penalty=a.mix_penalty,mix_normalization=a.mix_normalization,effective_mix_coefficient=mix_coefficient)),flush=True)
for step in range(a.steps):
 def pen(g,mu):return (torch.relu(mu+rho*g).square()-mu*mu)/(2*rho)
 def closure():
  opt.zero_grad();e2,gap,lam,score=vals();ge=e2/(a.cap*a.cap)-1;gg=(3.9155-gap)/3.9155
  loss=score/3000+pen(ge,mue)+pen(gg,mug)
  if a.mix_penalty:loss=loss+mix_coefficient*mix_metrics()[0]
  loss.backward();return loss
 if a.optimizer=='lbfgs':opt.step(closure)
 else:closure();opt.step()
 if step%(1 if a.optimizer=='lbfgs' else 20)==0:
  with torch.no_grad():ee,g,l,s=vals()
  er=float(ee.clamp(min=0).sqrt());sv=float(s.clamp(min=0).sqrt());gv=float(g)
  with torch.no_grad():mx,l1sq,trsq=mix_metrics()
  if er<=.019999 and gv>=3.9151 and sv<best:
   best=sv;save_fit(a.out,er)
  if er<=.019999 and gv>=3.9151 and float(mx)<bestmix:
   bestmix=float(mx)
   if a.mix_penalty:save_fit(str(Path(a.out).with_suffix(''))+'_lowmix.npz',er)
  mue=max(0,mue+rho*(float(ee)/(a.cap*a.cap)-1));mug=max(0,mug+rho*(3.9155-gv)/3.9155)
  if step%(5 if a.optimizer=='lbfgs' else 100)==0:print(json.dumps(dict(step=step,error_2=er,gap=gv,Lambda=float(l),lambda_eff=sv,best=best,
                                     mue=mue,mug=mug,mixing=float(mx),l1_squared=float(l1sq),trace_squared=float(trsq),bestmix=bestmix,mix_penalty=a.mix_penalty,seconds=time.time()-start)),flush=True)
 if a.optimizer=='adam':
  factor=.1+.9*(1-step/a.steps);opt.param_groups[0]['lr']=a.lr*factor;opt.param_groups[1]['lr']=a.offset_lr_factor*a.lr*factor
