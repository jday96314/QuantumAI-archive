"""Exact artifact, word, native-digest and immutable-reference audit."""
from core import *
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.scoreinfra.circuits import submission_sha256,load_circuits

def main():
    checked=0;vectors=0;maxdecode=0.;maxprojector=0.;graphs=set();binary=set()
    for folder in ['candidates','trajectories']:
        for path in sorted((HERE/folder).rglob('*.npz')):
            z=np.load(path)
            if 'angle_words' not in z or 'tree_edges' not in z:continue
            tr=z['tree_edges'];validate(tr);bits=int(z.get('brot',8));words=z['angle_words'].reshape(-1,75)
            assert np.issubdtype(words.dtype,np.integer) and np.all((words>=0)&(words<1<<bits))
            own=decode(words*np.pi/(1<<(bits-1)),tr);host=reconstruct_u_vectors(words/(1<<(bits-1)),tr)
            difference=float(np.max(abs(own-host)));assert difference<1e-12
            maxdecode=max(maxdecode,difference);checked+=1;vectors+=len(words);graphs.add(fingerprint(tr));binary.add(binary_signature(tr))
            if 'U' in z:
                U=z['U'];U=U/np.linalg.norm(U,axis=-1,keepdims=True)
                dot=np.abs((U.reshape(-1,N)*host[-U.shape[0]*U.shape[1]:]).sum(-1))
                error=float(np.max(abs(dot-1)));assert error<1e-11;maxprojector=max(maxprojector,error)
    natives=[]
    for file,source in [('incumbent_native.json',INC.parent),
        ('grid9_b351_submission_native.json',HERE/'candidates/grid9_b351_submission'),
        ('repair_topology_b351_submission_native.json',HERE/'candidates/repair_topology_b351_submission'),
        ('repair_nni_b479_submission_native.json',HERE/'candidates/repair_nni_b479_submission')]:
        r=json.loads((HERE/'results'/file).read_text());src=source/'submission.npz';tr=np.load(src)['tree_edges']
        digest=submission_sha256(load_circuits(source),src,tr);assert digest==r['submission_sha256']
        resource_match=None
        if source!=INC.parent:
            physical=json.loads((HERE/f'results/{source.name.removesuffix("_submission")}_physical.json').read_text())
            for key in ['toffoli','qubits','angles_qubits']:assert physical['resources'][key]==r['resources'][key]
            for key in ['error_1','error_2','Lambda','E_gap','lambda_eff']:assert abs(physical['metrics'][key]-r['metrics'][key])<1e-8
            resource_match=True
        natives.append(dict(record=file,status=r['status'],digest=digest,matches=True,verified_shots=r['verified_shots'],native_physical_resource_match=resource_match))
    baseline=json.loads((HERE/'results/incumbent_native.json').read_text());assert baseline['score']==708547388.9553173
    assert baseline['submission_sha256']=='ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7'
    sources={name:dict(path=str(path),sha256=sha(path)) for name,path in SOURCES.items()}
    hosts=['verify.py','scoreinfra/dfthc_math.py','scoreinfra/circuits.py','scoreinfra/costs_prepare.py','scoreinfra/reference.py']
    for relative in hosts:
        path=ROOT/'from_host/verifiers/dfthcsos'/relative;sources['host/'+relative]=dict(path=str(path),sha256=sha(path))
    report(HERE/'results/reference_provenance.json',sources,False)
    report(HERE/'results/artifact_audit.json',dict(checked_checkpoints=checked,decoded_vectors=vectors,
        distinct_tree_graphs=len(graphs),distinct_binary_lattices=len(binary),
        maximum_independent_decoder_difference=maxdecode,maximum_stored_projector_error=maxprojector,
        quarantined_invalid_encoder_outputs_excluded=True,native_digests=natives,incumbent_preserved=True))
    physical=[]
    for path in sorted((HERE/'results').glob('*_physical.json')):
        r=json.loads(path.read_text());physical.append(dict(candidate=path.stem,shape=r['shape'],bits=r['bits'],
             metrics=r['metrics'],physical_pass=r['physical_pass'],resources=r.get('resources'),conditional_product=r.get('conditional_product')))
    report(HERE/'results/physical_summary.json',physical,False)
    report(HERE/'selection.json',dict(replacement_selected=False,incumbent=str(INC.parent),
        incumbent_cost=baseline['score'],incumbent_native_record='results/incumbent_native.json',
        target_1_20_cost=baseline['score']/1.2,target_1_30_cost=baseline['score']/1.3,
        reason='Every new frozen submission fails two-body accuracy; measured normalization also exceeds cost ceilings.',
        classification=dict(A='No verified improvement',B='Validated encoders, lattice equivalence and modest conditional topology gains; no breakthrough claim',
            C='Tested topology/phase/repair pipelines uncompetitive, with fixed-dictionary certificates',
            D='Different well-conditioned exact-word dictionaries and broader shared-tree searches remain unresolved')))
    import scipy,platform
    report(HERE/'results/runtime_versions.json',dict(python=platform.python_version(),numpy=np.__version__,
        scipy=scipy.__version__,torch=torch.__version__,GPU=torch.cuda.get_device_name(0),
        torch_cpu_threads=torch.get_num_threads(),float_dtype='float64',tf32=False,
        timing_limitations='Independent jobs overlapped; wall times are descriptive, proposal/refit-work comparisons are authoritative'),False)
if __name__=='__main__':main()
