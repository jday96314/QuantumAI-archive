"""Translate and verify the syndrome-consuming BM circuit."""
import argparse,json,subprocess
from pathlib import Path
import kickmix as km
import published_dqi as p
from tower_inverse import apply_map,map_ops

def build(m,output,phase=False,sigma_test=False,tape_prefix='stream'):
    mod=p.load();forward,back,transform=p.isomorphism(mod)
    hook=json.loads((p.HERE/f'parameters/{tape_prefix}_hook.json').read_text())
    b=km.CircuitBuilder();syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    masks=[b.create_classical_register(12,f'E{i}') for i in range(m)]
    q=[x for w in syn for x in w]+list(b.alloc_qubits(hook['maxq']-840))
    for w in syn:apply_map(b,w,map_ops(forward))
    b.info_tag('SPLIT_START')
    if sigma_test:
        for w in hook['sigma']:b.z(km.array([q[i] for i in w]))
    if phase:
        from stream_forney import phase_kernel
        b.free(km.array([q[i] for i in hook['free']]))
        phase_kernel(b,hook,masks,q,forward,back)
    b.info_tag('SPLIT_END')
    for w in syn:apply_map(b,w,map_ops(forward),True)
    c=b.finish_circuit();wrapper=output.with_suffix('.wrapper.kmx');c.to_file(str(wrapper),format='kmx')
    paths=[output.with_suffix('.'+part+'.kmx') for part in ['pre','sweep','post']]
    handles=[x.open('w') for x in paths];part=0
    for line in wrapper.open():
        if 'SPLIT_START' in line:part=1
        elif 'SPLIT_END' in line:part=2
        else:handles[part].write(line)
    for f in handles:f.close()
    full=output.with_suffix('.kmx')
    command=[str(p.HERE/'translate_published'),str(p.HERE/f'parameters/{tape_prefix}_tape.npy'),str(hook['index']),
                    str(c.num_bits),*map(str,paths),str(full),str(hook['outcomes'])]
    if 'phases' in hook:
        phases=output.with_suffix('.phases.txt')
        phases.write_text('\n'.join(str(i)+' '+name for i,name in hook['phases'])+'\n')
        command.append(str(phases))
    subprocess.run(command,check=True)
    c=km.Circuit.from_file(str(full));print('built',c.num_qubits,c.num_bits,c.operation_statistics().num_toffoli_like,flush=True)
    c.to_file(str(output),format='kmb')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--m',type=int,default=71);a.add_argument('--phase',action='store_true')
    a.add_argument('--sigma-test',action='store_true')
    a.add_argument('--tape-prefix',default='stream')
    a.add_argument('--output',type=Path,required=True);args=a.parse_args();build(args.m,args.output,args.phase,args.sigma_test,args.tape_prefix)
