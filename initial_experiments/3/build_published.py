"""Fast, auditable IR translation through kickmix's native text parser."""
import argparse,json,subprocess,time
from pathlib import Path
import kickmix as km
import published_dqi as p
from tower_inverse import map_ops,apply_map

def build(m,output):
    pub=p.load();forward,back,transform=p.isomorphism(pub)
    p.ISOM_R=pub.LOG[forward[1]]
    hook=json.loads((p.HERE/'parameters/published_hook.json').read_text())
    b=km.CircuitBuilder()
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    masks=[b.create_classical_register(12,f'E{i}') for i in range(m)]
    q=[v for w in syn for v in w]+list(b.alloc_qubits(hook['maxq']-840))
    ops=map_ops(forward)
    for w in syn:apply_map(b,w,ops)
    b.info_tag('SPLIT_SWEEP_START')
    p.sweep(b,pub,hook,masks,q,back,transform,range(m))
    b.info_tag('SPLIT_SWEEP_END')
    for w in syn:apply_map(b,w,ops,True)
    c=b.finish_circuit();wrapper=output.with_suffix('.wrapper.kmx')
    c.to_file(str(wrapper),format='kmx')
    paths=[output.with_suffix('.'+x+'.kmx') for x in ['pre','sweep','post']]
    handles=[x.open('w') for x in paths];part=0
    with wrapper.open() as f:
        for line in f:
            if 'SPLIT_SWEEP_START' in line:part=1
            elif 'SPLIT_SWEEP_END' in line:part=2
            else:handles[part].write(line)
    for f in handles:f.close()
    full=output.with_suffix('.kmx')
    subprocess.run([str(p.HERE/'translate_published'),str(p.HERE/'parameters/published_tape.npy'),
                    str(hook['index']),str(c.num_bits),*map(str,paths),str(full),str(hook['outcomes'])],check=True)
    c=km.Circuit.from_file(str(full))
    print('built',c.num_qubits,c.num_bits,c.operation_statistics().num_toffoli_like,flush=True)
    c.to_file(str(output),format='kmb')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--m',type=int,default=4095);a.add_argument('--output',type=Path,required=True)
    args=a.parse_args();build(args.m,args.output)
