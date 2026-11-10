"""Low et al. residual + Lambda + gap objective at the native fixed shift.

From-scratch random unit rows/Gaussian weights; no incumbent float basin is
used. Optional expected error over independent uniform half-LSB angle noise
uses ordinary reparameterized gradients, never straight-through rounding.
All published checkpoints recompute the exact full native signed metric.
"""
import argparse,time,json,copy,hashlib,os
from audit_core import *
HERE=Path(__file__).resolve().parent

def continuous_encode(u,tree):
    v=u.copy().reshape(-1,76);out=np.zeros((len(v),75))
    for j in reversed(range(75)):
        p,c=tree[j];a,b=v[:,p].copy(),v[:,c].copy();sg=np.where(a<0,-1.,1.)
        out[:,j]=np.arctan2(b*sg,a*sg)*128/np.pi
        v[:,p]=sg*np.hypot(a,b);v[:,c]=0
    return out.reshape(*u.shape[:-1],75)

def main():
    p=argparse.ArgumentParser();p.add_argument('--R',type=int,default=1);p.add_argument('--B',type=int,default=719);p.add_argument('--C',type=int,default=350)
    p.add_argument('--seed',type=int,default=1);p.add_argument('--steps',type=int,default=20000)
    p.add_argument('--lambda-reg',type=float,default=300);p.add_argument('--epsilon-reg',type=float,default=.02)
    p.add_argument('--lr',type=float,default=.01);p.add_argument('--final-lr',type=float,default=1e-5)
    p.add_argument('--tag',required=True);p.add_argument('--checkpoint',type=int,default=2000)
    p.add_argument('--init');p.add_argument('--noise',type=int,default=0,help='number of half-LSB noise samples per expected-error step')
    p.add_argument('--weight-cap',type=float,default=100.);p.add_argument('--lambda-cap',type=float,default=320.)
    p.add_argument('--hard-lambda-penalty',type=float,default=0.)
    p.add_argument('--precision',choices=['float64','float32'],default='float64',help='float32 uses a stable full residual; checkpoints and physical quantities remain float64')
    p.add_argument('--enforce-lambda',action='store_true',help='reject/backtrack updates crossing the Lambda cap; requires a feasible initializer')
    p.add_argument('--checkpoint-device',choices=['cpu','same'],default='cpu')
    a=p.parse_args();torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed)
    source_hashes={}
    for source in [Path(__file__).resolve(),HERE.parent.parent/'5/chem-audit/audit_core.py']:
        raw=source.read_bytes();digest=hashlib.sha256(raw).hexdigest();source_hashes[source.name]=digest
        archive=HERE/'results/code';archive.mkdir(exist_ok=True)
        (archive/f'{digest}_{source.name}').write_bytes(raw)
    base=dict(np.load(INC/'submission.npz'));tree=base['tree_edges']
    train_device=os.environ.get('AUDIT_DEVICE','cuda')
    m=Metric(tree,device='cpu' if a.checkpoint_device=='cpu' else train_device);R,B,C=a.R,a.B,a.C
    if a.init:
        z,u,w,tree,words=load(a.init);R,B,C=w.shape
        off=z['W'] if z['W'].ndim==2 else z['H2_I']-w.sum(1)
    else:
        u=rng.normal(size=(R,B,76));u/=np.linalg.norm(u,axis=-1,keepdims=True)
        w=rng.normal(scale=.01,size=(R,B,C));off=-w.sum(1)
    mt=copy.copy(m);mt.d=dict(device=train_device,dtype=torch.float32 if a.precision=='float32' else torch.float64)
    mt.G=m.G.to(**mt.d);mt.fac=m.fac.to(**mt.d)
    mt.ii=m.ii.to(train_device);mt.jj=m.jj.to(train_device)
    mt.stages=[tuple(x.to(train_device) for x in stage) for stage in m.stages]
    mt.h0=m.h0.to(device=train_device,dtype=torch.float64)
    coord=mt.tensor(continuous_encode(u,tree) if a.noise else u).requires_grad_()
    wt=mt.tensor(w).requires_grad_();wb=mt.tensor(off+w.sum(1)).requires_grad_()
    with torch.no_grad():wt.mul_(min(1.,a.weight_cap/max(float(wt.norm()),1e-20)))
    # Angles are measured in word units, requiring a larger coordinate step.
    coord_factor=20. if a.noise else 1.
    opt=torch.optim.Adam([{'params':[coord],'lr':coord_factor*a.lr},{'params':[wt,wb],'lr':a.lr}])
    const=float(m.G.square().sum());start=time.time();rows=[];best=float('inf');best_valid=float('inf');cap_backtracks=0
    best_integer_error=float('inf');best_integer_physical=float('inf');best_round_physical=float('inf')
    def unpack_state():
        uu=mt.decode(coord) if a.noise else coord/coord.norm(dim=-1,keepdim=True)
        return uu,wb.double()-wt.double().sum(1)
    def error2(uu):
        P=mt.feats(uu);X=(wt.transpose(1,2)@P).flatten(0,1)
        if a.precision=='float32':return (X.T@X-mt.G).square().sum().clamp(min=1e-20)
        # Algebraically the full signed native metric; no PSD clipping.
        return (const-2*(X*(X@mt.G)).sum()+(X@X.T).square().sum()).clamp(min=1e-20)
    def metrics(uu,oo):
        un=uu.double();un=un/un.norm(dim=-1,keepdim=True);ww=wt.double();oo=oo.double()
        coefficients=(ww@oo[...,None]).squeeze(-1)
        hc=torch.bmm(un.transpose(1,2),un*coefficients[...,None]).sum(0)
        eig=torch.linalg.eigvalsh(mt.h0-hc);ab=eig.abs().sum()
        gap=m.constant+ab-eig.sum()+.5*oo.square().sum()
        lam=ab+.25*(ww.abs().sum(1)+wb.double().abs()).square().sum()
        return gap,lam
    for step in range(a.steps+1):
        if step%a.checkpoint==0 or step==a.steps:
            with torch.no_grad():
                uu,oo=unpack_state();cp,quant=checkpoint(m,uu.to(**m.d),wt.to(**m.d),oo.to(**m.d),tree,2)
                train_squared=float(error2(uu));contraction_discrepancy=abs(train_squared-cp['float']['error_2']**2)
                norm_discrepancy=abs(np.sqrt(train_squared)-cp['float']['error_2'])
                # Check norm units: an absolute squared-error tolerance becomes
                # disproportionately strict as the residual shrinks.
                if a.precision=='float32':assert norm_discrepancy<1e-5,norm_discrepancy
                else:assert contraction_discrepancy<1e-8,contraction_discrepancy
                g,l=metrics(uu,oo)
                physical_discrepancy=max(abs(float(g)-cp['float']['E_gap']),abs(float(l)-cp['float']['Lambda']))
                assert physical_discrepancy<1e-7,physical_discrepancy
                objective=cp['float']['error_2']/a.epsilon_reg+float(l)/a.lambda_reg+max(0.,(float(g)-3.91515)/3.91515)
                rec=dict(step=step,shape=[R,B,C],objective=objective,metrics=cp,seconds=time.time()-start,
                    squared_error_contraction_discrepancy=contraction_discrepancy,cap_backtracks=cap_backtracks)
                rec['training_norm_discrepancy']=norm_discrepancy
                rec['physical_contraction_discrepancy']=physical_discrepancy
                if a.enforce_lambda:assert cp['float']['Lambda']<=a.lambda_cap+1e-6
                if a.noise:
                    samples=[]
                    for _ in range(8):
                        un=mt.decode(coord+torch.rand_like(coord)-.5)
                        gn,ln=metrics(un,oo)
                        samples.append(dict(error_2=float(m.loss(un.to(**m.d),wt.to(**m.d)).sqrt()),Lambda=float(ln),E_gap=float(gn)))
                    rec['noise_samples']=samples
                    rec['expected_error_2']=float(np.mean([x['error_2'] for x in samples]))
                    rec['float_objective']=objective
                    objective=rec['expected_error_2']/a.epsilon_reg+float(l)/a.lambda_reg+max(0.,(float(g)-3.91515)/3.91515)
                    rec['objective']=objective
                rows.append(rec);print(json.dumps(rec),flush=True)
                out=HERE/f'parameters/{a.tag}.npz'
                if objective<best:
                    best=objective;save(out,uu,wt,tree,W=oo.cpu().numpy(),err2=cp['float']['error_2'])
                if cp['float']['error_2']<.02 and cp['float']['Lambda']<best_valid:
                    best_valid=cp['float']['Lambda'];save(HERE/f'parameters/{a.tag}_passing_float.npz',uu,wt,tree,W=oo.cpu().numpy())
                uq,wq,oq,aw=quant
                if cp['refit']['error_2']<best_integer_error:
                    best_integer_error=cp['refit']['error_2']
                    save(HERE/f'parameters/{a.tag}_best_integer.npz',uq,wq,tree,aw,W=oq.cpu().numpy())
                for kind,ww,oo2 in [('refit',wq,oq),('rounded',wt.double(),oo)]:
                    row=cp[kind]
                    if row['error_2']<=.02 and row['E_gap']>3.914898 and row['Lambda']>49 and row['lambda_eff'] is not None:
                        old=best_integer_physical if kind=='refit' else best_round_physical
                        if row['lambda_eff']<old:
                            if kind=='refit':best_integer_physical=row['lambda_eff']
                            else:best_round_physical=row['lambda_eff']
                            save(HERE/f'parameters/{a.tag}_best_physical_{kind}.npz',uq,ww,tree,aw,W=oo2.detach().cpu().numpy())
                save(HERE/f'parameters/{a.tag}_rounded.npz',uq,wt.double(),tree,aw,W=oo.cpu().numpy(),err2=cp['rounded']['error_2'])
                save(HERE/f'parameters/{a.tag}_integer.npz',uq,wq,tree,aw,W=oq.cpu().numpy(),err2=cp['refit']['error_2'])
                save(HERE/f'parameters/{a.tag}_last.npz',uu,wt,tree,W=oo.cpu().numpy())
                write_json(HERE/f'results/{a.tag}.json',dict(arguments=vars(a),training_device=train_device,checkpoint_device=m.d['device'],
                    source_sha256=source_hashes,rows=rows,finished=step==a.steps))
        if step==a.steps:break
        opt.zero_grad();uu,oo=unpack_state();g,l=metrics(uu,oo)
        if a.noise:
            for sample in range(a.noise):
                perturbation=torch.rand_like(coord)-.5 if sample%2==0 else -perturbation
                # Linearity of differentiation permits streaming the expected
                # residual gradient without retaining all noisy N^4 graphs.
                (error2(mt.decode(coord+perturbation)).sqrt()/(a.noise*a.epsilon_reg)).backward()
            loss=l/a.lambda_reg+torch.relu((g-3.91515)/3.91515)
        else:loss=error2(uu).sqrt()/a.epsilon_reg+l/a.lambda_reg+torch.relu((g-3.91515)/3.91515)
        loss=loss+a.hard_lambda_penalty*torch.relu(l/a.lambda_cap-1).square()
        loss.backward()
        if a.enforce_lambda:
            old=[v.detach().clone() for v in [coord,wt,wb]]
        opt.step()
        with torch.no_grad():
            wt.mul_((a.weight_cap/wt.norm().clamp(min=1e-20)).clamp(max=1.))
            if not a.noise:coord.div_(coord.norm(dim=-1,keepdim=True))
            if a.enforce_lambda:
                un,on=unpack_state()
                if float(metrics(un,on)[1])>a.lambda_cap:
                    new=[v.detach().clone() for v in [coord,wt,wb]];feasible=False
                    for scale in [.5,.25,.125,.0625,.015625,.00390625]:
                        cap_backtracks+=1
                        for v,previous,proposal in zip([coord,wt,wb],old,new):v.copy_(previous+scale*(proposal-previous))
                        un,on=unpack_state()
                        if float(metrics(un,on)[1])<=a.lambda_cap:feasible=True;break
                    if not feasible:
                        for v,previous in zip([coord,wt,wb],old):v.copy_(previous)
        lr=a.lr*(a.final_lr/a.lr)**((step+1)/a.steps)
        opt.param_groups[0]['lr']=coord_factor*lr;opt.param_groups[1]['lr']=lr

if __name__=='__main__':main()
