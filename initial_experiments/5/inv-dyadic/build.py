"""Independent adaptive-prefix inverse-square-root circuit.

Both architectures use integer tables from fit.cpp. 'local' recenters input
bits reversibly at variable depths; 'shared' retains one coarse coordinate.
The only reused quantum components are normalization, radix-four multiply-add,
the unequal-width adder, and measured ANF erasure. No incumbent builder runs.
"""
import argparse, collections, hashlib, json, subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt import nr_reference_circuit as nr
from binary_centered_invsqrt import cshift as normal_shift, del_flag_if_zero
from short_adder import add
from booth_adder import accumulate
from monomial_unlookup import erase
from certify import certificate
from recenter import transform

HERE=Path(__file__).resolve().parent
OLD=HERE.parents[1]/'4/codex'

class ProfileBuilder:
    def __init__(self):
        self.inner=km.CircuitBuilder();self.stage='start';self.events=[];self.peaks=collections.defaultdict(int)
    def __getattr__(self,name):return getattr(self.inner,name)
    def sample(self):
        n=self.num_allocated_qubits;self.peaks[self.stage]=max(self.peaks[self.stage],n)
    def alloc_qubits(self,n):
        r=self.inner.alloc_qubits(n);self.sample();return r
    def create_quantum_register(self,n,name):
        r=self.inner.create_quantum_register(n,name);self.sample();return r
    def free(self,r):self.inner.free(r)
    def mark_stage(self,name):
        self.stage=name;self.sample();self.events.append(dict(stage=name,live=self.num_allocated_qubits));self.info_tag(name)

def tree_xor(b,address,targets,rows):
    """A compressed prefix tree: one AND per nonconstant non-root node."""
    rows=np.asarray(rows,dtype=np.int64)
    def emit(ctrl,row):
        for target,v in zip(targets,row):
            if len(target):b.gf2_iadd(int(v),target=target,control=ctrl)
    def visit(ctrl,level,lo,hi):
        if np.all(rows[lo:hi]==rows[lo]):emit(ctrl,rows[lo]);return
        assert level>=0
        bit=address[level];mid=(lo+hi)//2;t=b.alloc_qubits(1)[0]
        b.x(bit);b.init_and(ctrl,bit,t);b.x(bit)
        visit(t,level-1,lo,mid);b.cx(ctrl,t);visit(t,level-1,mid,hi)
        b.del_and(ctrl,bit,t);b.free(km.array([t]))
    mid=len(rows)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid)
    b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(rows))

def build(model,limit=90,keep=4,delay=2,product_limit=17,deferred=3,split_delta=0,stream_erase=False):
    b=ProfileBuilder();wc,w1,w2=model['widths']
    local=model['mode'].startswith('local');base=min(l[2] for l in model['leaves'])
    depth=max(l[2] for l in model['leaves']);K=depth-base if local else 0
    # Shared-coordinate fits may start at a depth below their actual leaves.
    coord_base=min(l[3] for l in model['leaves'])
    table=np.zeros((1<<(depth+1),4),dtype=np.int64)
    for p,lo,d,f,degree,a,slope,curve,error in model['leaves']:
        q=lo>>(21-base)
        residual=a+(q<<(wc-base-2))-(11<<(wc-4))-(p<<(wc-2))
        begin=(p<<depth)+(lo>>(21-depth));end=begin+(1<<(depth-d))
        table[begin:end]=[residual,slope,curve if degree==2 else 0,(1<<(d-base))-1 if local else 0]
    residual_bias=0
    if table[:,0].min()<0:
        residual_bias=1<<(wc-7)
        while table[:,0].min()+residual_bias<0:residual_bias*=2
        table[:,0]+=residual_bias
    sizes=[int(table[:,j].max()).bit_length() for j in range(3)]
    assert sizes[0]<=wc-3
    bound=model.get('affine_residual_range')
    integrate_full=bound is None or bound[0]+residual_bias<0 or bound[1]+residual_bias>=(1<<sizes[0])
    split=[wc-base-2+split_delta,w1-base-2+split_delta]
    assert all(0<sp<sz for sp,sz in zip(split,sizes))
    x=b.create_quantum_register(22,'x');flags=b.alloc_qubits(5)
    def normalize(inv):
        for k in (range(5) if inv else reversed(range(5))):
            shift=1<<k
            if inv:
                normal_shift(b,x,shift,flags[k],22-shift,inv=True);del_flag_if_zero(b,x[-shift:],flags[k])
            else:
                b.flip_if_equal(x[-shift:],0,target=flags[k]);normal_shift(b,x,shift,flags[k],22-shift)
    b.mark_stage('normalize');normalize(False)
    saved=b.alloc_qubits(K)
    for j,q in enumerate(saved):b.cx(x[20-base-j],q)
    address=km.array([*x[21-depth:21],flags[0]])
    if K:
        address=km.array([*(saved[depth-1-base-i] if base<=depth-1-i<depth else x[21-depth+i] for i in range(depth)),flags[0]])
    depth_flags=b.alloc_qubits(K)
    core=b.alloc_qubits(sizes[0]);slope=b.alloc_qubits(sizes[1]);curve=b.alloc_qubits(sizes[2])
    regs=[core,slope,curve,depth_flags]
    cutoff=sizes[0]-delay;assert cutoff>=split[0]
    lowcore=core[:cutoff];delayed=core[cutoff:]
    loadrows=table.copy();loadrows[:,0]&=(1<<cutoff)-1
    b.mark_stage('lookup');tree_xor(b,address,[lowcore,slope,curve,depth_flags],loadrows)

    def recenter(inv):
        if not local:
            b.x(x[20-coord_base]);return
        transform(b,x,saved,depth_flags,base,inv)
    b.mark_stage('coordinate');recenter(False)
    # Private phase headroom grows only after high curve bits have been erased.
    initial_extra=wc+9-sizes[0]-deferred
    extra=b.alloc_qubits(initial_extra)
    private=km.array([*core,*extra]);y=km.array([*private[:wc-1],x[21],*private[wc-1:]])
    counters=[b.alloc_qubits(1),b.alloc_qubits(1)]
    def workspace():return km.array([q for q in y[sizes[0]:] if q!=x[21]])
    delayed_measure=None
    def integrate(j,inv):
        reg=(y[:wc-1] if integrate_full else core) if j==0 else slope;high=reg[split[j]:]
        sign=x[20-coord_base];b.x(sign);b.cx(sign,high);b.x(sign)
        carry_work=y[wc:] if j==0 else workspace()
        add(b,counters[j],high,subtract=inv,carry_workspace=carry_work,qubit_limit=limit)
        b.x(sign);b.cx(sign,high);b.x(sign)
    def after_coordinate_lookup(targets,rows):
        flip=not local and depth>coord_base
        if flip:b.x(x[20-coord_base])
        tree_xor(b,address,targets,rows)
        if flip:b.x(x[20-coord_base])
    def horner(j,inv):
        nonlocal delayed_measure
        target=core if j==0 else slope;source=slope if j==0 else curve
        target=km.array([*target[:split[j]],*counters[j]])
        if inv:
            integrate(j,True)
            if j==0 and delay:
                delayed_measure=b.alloc_bits(delay);b.hmr(delayed,delayed_measure)
        b.x(x[21])
        accumulate(b,source,x[:21-coord_base],target,(wc-w1 if j==0 else w1-w2)-21,
                   subtract=not inv,rounded=True,
                   borrow=km.array([*workspace(),*delayed]),auxiliary=km.array([x[21]]),
                   qubit_limit=limit,product_limit=product_limit)
        b.x(x[21])
        if not inv:
            if j==0 and delay:
                # Lookup address is preserved independently of recentered bits.
                after_coordinate_lookup([delayed],((table[:,0]>>cutoff)&((1<<delay)-1))[:,None])
            integrate(j,False)
    b.mark_stage('curve_horner');horner(1,False)
    keep=min(keep,len(curve));high_curve=curve[keep:];curve_measure=b.alloc_bits(len(high_curve))
    b.mark_stage('release_curve');b.hmr(high_curve,curve_measure);b.free(high_curve)
    later=b.alloc_qubits(deferred);private=km.array([*private,*later]);y=km.array([*y,*later])
    b.mark_stage('linear_horner');horner(0,False)

    def restore(inv):
        head=y[wc:]
        def constant(inv):
            a,c,d=y[wc-4:wc-1];scratch=head[0]
            if not inv:b.x(a);b.cx(a,c)
            b.init_and(a,c,scratch);b.x(scratch);b.cx(scratch,d);b.x(scratch);b.del_and(a,c,scratch)
            if inv:b.cx(a,c);b.x(a)
        def parity(inv):b.cx(flags[0],y[wc-2])
        def bias(inv):
            if not residual_bias:return
            bit=residual_bias.bit_length()-1;one=head[-1];b.x(one)
            add(b,km.array([one]),y[bit:wc-1],subtract=not inv,carry_workspace=head[:-1],qubit_limit=limit)
            b.x(one)
        def coarse(inv):add(b,x[21-base:21],y[wc-base-2:wc-1],subtract=not inv,carry_workspace=head,qubit_limit=limit)
        # The original coarse address may overlap the coordinate in local mode.
        def coarse_original(inv):
            if not K:return coarse(inv)
            q=km.array([*x[21-base:21]])
            add(b,q,y[wc-base-2:wc-1],subtract=not inv,carry_workspace=head,qubit_limit=limit)
        steps=[constant,parity,coarse_original,bias]
        for step in (reversed(steps) if inv else steps):step(inv)
    def shift(inv):
        occ=wc
        stages=[]
        for k in range(4):stages.append((1<<k,occ,flags[k+1]));occ+=1<<k
        for step,occupied,ctrl in (reversed(stages) if inv else stages):
            indices=list(range(min(occupied,len(y)-step)))
            for i in (indices if inv else reversed(indices)):
                src,dst=y[i],y[i+step]
                if i==wc-1 and step==1:
                    b.cx(ctrl,src);b.cx(ctrl,dst)
                elif i+step>=occupied:
                    if inv:b.cx(dst,src);b.del_and(ctrl,src,dst)
                    else:b.init_and(ctrl,src,dst);b.cx(dst,src)
                else:b.cswap(ctrl,src,dst)
    b.mark_stage('restore_affine');restore(False)
    b.mark_stage('phase');shift(False);nr._phase_by_bits(b,y);shift(True)
    b.mark_stage('remove_affine');restore(True)
    b.mark_stage('inverse_linear');horner(0,True)
    b.free(later);private=private[:-deferred] if deferred else private;y=y[:-deferred] if deferred else y
    high_curve=b.alloc_qubits(len(high_curve));curve=km.array([*curve[:keep],*high_curve]);regs[2]=curve
    b.mark_stage('replay_curve');after_coordinate_lookup([high_curve],(table[:,2]>>keep)[:,None])
    for bit,q in zip(curve_measure,high_curve):
        b.push_condition(bit);b.z(q);b.pop_condition()
    b.bit_store(curve_measure,False);b.free(curve_measure)
    b.mark_stage('inverse_curve');horner(1,True)
    b.free(extra)
    for counter in counters:b.free(counter)
    b.mark_stage('restore_coordinate');recenter(True)
    # Measurements of delayed c0 high bits join the final coefficient erasure.
    targets=[lowcore,slope,curve,depth_flags]
    widths=[len(r) for r in targets];packed=[]
    for row in loadrows:
        value=sum(int(v)<<sum(widths[:j]) for j,v in enumerate(row))
        if delay:value|=(int(table[len(packed),0])>>cutoff)<<sum(widths)
        packed.append(value)
    erase_fn=erase
    target=km.array([q for r in targets for q in r]);needed=sum((1<<n)-n-1 for n in [len(address)//2,len(address)-len(address)//2])
    if stream_erase:
        from stream_erase import erase as erase_fn,workspace_size
        needed=workspace_size(len(address))
    # Append spare measured zeros before delayed outcomes, shifting table bits.
    spare=b.alloc_qubits(max(0,needed-len(target)))
    if len(spare):
        packed=[(v&((1<<len(target))-1))|((v>>len(target))<<(len(target)+len(spare))) for v in packed]
        target=km.array([*target,*spare])
    b.mark_stage('erase_lookup');erase_fn(b,address,target,packed,extra_outcomes=delayed_measure if delay else ())
    b.free(spare)
    for r in regs:b.free(r)
    for j,q in enumerate(saved):b.cx(x[20-base-j],q)
    b.free(saved)
    b.mark_stage('denormalize');normalize(True);b.free(flags)
    return b.finish_circuit(),dict(stages=b.events,stage_live_peaks=dict(b.peaks),sizes=sizes,
             split=split,base=base,max_depth=depth,coordinate_base=coord_base,coordinate_copies=K,
             residual_bias=residual_bias,integrate_full=integrate_full,
             parameters=dict(limit=limit,keep=keep,delay=delay,product_limit=product_limit,deferred=deferred,split_delta=split_delta,stream_erase=stream_erase))

def save(c,stem,stage_names):
    raw=HERE/'candidates'/f'{stem}.raw.kmx';txt=raw.with_suffix('.compact.kmx');path=raw.with_name(stem+'.kmb')
    c.to_file(str(raw),format='kmx')
    run=subprocess.run([str(OLD/'compact_wires'),str(raw),str(txt),'22'],check=True,capture_output=True,text=True)
    c=km.Circuit.from_file(str(txt));c.to_file(str(path),format='kmb')
    # Keep tagged source for review and component accounting.
    components={};stage='start'
    for line in raw.read_text().splitlines():
        if line.startswith('INFO_TAG '):
            label=line.split('"')[1]
            if label in stage_names:stage=label
        elif line.split() and line.split()[0] in ['CCX','CCZ','CSWAP']:components[stage]=components.get(stage,0)+1
    raw.unlink();txt.unlink()
    return c,path,dict(compaction=run.stderr.strip(),component_static_toffolis=components)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('model',type=Path);p.add_argument('--limit',type=int,default=90)
    p.add_argument('--keep',type=int,default=4);p.add_argument('--delay',type=int,default=2)
    p.add_argument('--product',type=int,default=17);p.add_argument('--deferred',type=int,default=3)
    p.add_argument('--split-delta',type=int,default=0);p.add_argument('--stream-erase',action='store_true')
    a=p.parse_args();model=json.loads(a.model.read_text())
    cert_path=HERE/'results'/f'{a.model.stem}_certificate.json'
    cert=json.loads(cert_path.read_text()) if cert_path.exists() else certificate(model)
    model.update(affine_residual_range=cert.get('affine_residual_range'))
    stem=f'{a.model.stem}_q{a.limit}_k{a.keep}_d{a.delay}_p{a.product}_f{a.deferred}_s{a.split_delta}'
    if a.stream_erase:stem+='_stream'
    c,report=build(model,a.limit,a.keep,a.delay,a.product,a.deferred,a.split_delta,a.stream_erase)
    c,path,resources=save(c,stem,report['stage_live_peaks'])
    report.update(resources);report.update(status='BUILT',model=str(a.model),path=str(path),
        Q=c.num_qubits,T_static=c.operation_statistics().num_toffoli_like,
        circuit_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),certificate=cert)
    report['conservative_cost']=c.num_qubits*(report['T_static']+cert['synthesis_toffoli_equivalent'])
    (HERE/'results'/f'{stem}_build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
