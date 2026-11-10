"""Integrate the validated repaired H map into the unchanged accepted core.

All prior source files remain read-only. The local generated build source
differs only at the phase hook and the scan. Full native verification is
required before replacing the accepted artifact.
"""
import argparse,inspect,json,hashlib
from pathlib import Path
import kickmix as km
import dqi_budget_probe as probe
import repaired_sylvester as repaired
import low_workspace_inverse
from safe_retained_inverse import bind
from import_published import load
from published_dqi import isomorphism
from tower_inverse import map_ops,apply_map,linear_xor

HERE=Path(__file__).resolve().parent
h=repaired.fast.module('agent20_decoder_h_phase','h_phase.py')

def prepare_H(m,E,bp,bb,p,stash,phase_shift,flags):
    A=[*bp,p];B=[*bb,stash[1]]
    E.x(stash[1][0])
    repair=repaired.Repair(m,E,flags)
    repair.rotation(A,phase_shift,1)
    h.sylvester=repair.apply
    H=h.prepare(m,E,A,B)
    return dict(repair=repair,A=A,B=B,H=H,stash=stash,phase_shift=phase_shift)

def unprepare_H(m,E,state):
    repair=state['repair'];h.sylvester=repair.apply
    h.prepare(m,E,state['A'],state['B'],True)
    repair.rotation(state['A'],state['phase_shift'],1,True)
    repair.close();E.x(state['stash'][1][0])

SCAN='''    recycled=True
    H=[km.array([allq[q] for q in word]) for word in H_state['H']]
    B.append(km.array([allq[q] for q in stash[1]]))
    counter=km.array([allq[q] for q in H_state['repair'].count])
    def rotate_banks(back=False):
        for i in (reversed(range(len(counter))) if back else range(len(counter))):
            for bank in [H,B]:rotate(kb,bank,counter[i],1<<i,back)
    rotate_banks()
    root,value=B[0],H[0]
    gf=km.GF2Field(12)
    def lin(word,target,c):
        linear_xor(kb,word,target,[transform(m.gmul(c,1<<i),back) for i in range(12)])
    for j,mask in enumerate(masks):
        beta=gf.pow(2,(-j)%4095);x=transform(beta,forward)
        # The two accumulator coefficients are changed to host field basis;
        # all other polynomial words stay in the arithmetic emitter's basis.
        apply_map(kb,root,map_ops(back));apply_map(kb,value,map_ops(back))
        def evaluate():
            for i,word in enumerate(B[1:],1):lin(word,root,m.gpow(x,i))
            for i,word in enumerate(H[1:],1):lin(word,value,m.gpow(x,i))
        evaluate()
        constant=gf.pow(beta,N)
        columns=[sum(((gf.mul(constant,gf.mul(1<<z,1<<z))>>i)&1)<<z for z in range(12)) for i in range(12)]
        maskops=map_ops(columns);apply_map(kb,mask,maskops)
        h_phase(kb,root,value,mask)
        apply_map(kb,mask,maskops,True);evaluate()
        apply_map(kb,value,map_ops(back),True);apply_map(kb,root,map_ops(back),True)
        if j%256==0:print(json.dumps(dict(stage='H_scan',position=j,seconds=time.time()-begin)),flush=True)
    rotate_banks(True)
    kb.info_tag('COMPACT_SCAN_END')'''

def builder(profile=False,small=False):
    source=inspect.getsource(probe.decoder.build)
    old='    scale();hook=len(b._items)\n    scale(True)'
    new="    scale()\n    H_state=prepare_H(m,E,bp,bb,p,stash,phase_shift,flags)\n    H_held=[*H_state['repair'].tags,*H_state['repair'].gamma,*H_state['repair'].count]\n    hook=len(b._items)\n    unprepare_H(m,E,H_state)\n    scale(True)"
    assert source.count(old)==1;source=source.replace(old,new)
    old='    flagset=set(flags+waits+[terminal,phase_shift]+[q for w in stash for q in w])'
    assert source.count(old)==1;source=source.replace(old,old[:-1]+'+H_held)')
    start=source.index('    recycled=N//2>=35')
    end=source.index("    kb.info_tag('COMPACT_SCAN_END')",start)+len("    kb.info_tag('COMPACT_SCAN_END')")
    source=source[:start]+SCAN+source[end:]
    namespace=dict(probe.decoder.__dict__)
    namespace.update(prepare_H=prepare_H,unprepare_H=unprepare_H,h_phase=h.phase,
                     load=probe.count_module if profile else load,
                     BoundedArithmetic=probe.trimmed.TrimmedArithmetic,
                     frames=probe.frames)
    if small:
        from blocked_euclid_sequence import frames
        namespace['frames']=frames
    exec(source,namespace);build=namespace['build'];build.profile_only=profile
    build.slot_counts=None if small else probe.counts
    (HERE/'parameters/H_decoder_build_source.py').write_text(source)
    low_workspace_inverse.bind=bind
    return build

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--profile',action='store_true');p.add_argument('--small',action='store_true');a=p.parse_args()
    build=builder(a.profile,a.small);N,M,width=(14,127,2) if a.small else (70,4095,3)
    path=HERE/'candidates'/('repaired_H_small_n14.kmb' if a.small else 'repaired_H_dqi.kmb')
    report=build(N,width,100,M,path,linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
    if a.profile:
        (HERE/'results/repaired_H_core_profile.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
    else:
        manifest=dict(status='BUILT_REQUIRES_NATIVE_VERIFICATION',circuit=str(path),
                      sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      Q=report.num_qubits,static_T=report.operation_statistics().num_toffoli_like,
                      independent_salt='agent20-repaired-H-dqi-independent-20261025',
                      source_sha256={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),HERE/'repaired_sylvester.py',HERE/'fast_sylvester.py',HERE/'parameters/H_decoder_build_source.py']})
        (HERE/'results'/('repaired_H_small_build.json' if a.small else 'repaired_H_build.json')).write_text(json.dumps(manifest,indent=2)+'\n')
        print(json.dumps(manifest),flush=True)
