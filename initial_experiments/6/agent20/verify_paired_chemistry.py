"""All-address dirty-target checks, then two unchanged full host CLI runs."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import kickmix as km
import numpy as np
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.circuits import load_circuits

HERE = Path(__file__).resolve().parent
DIRECTORY = HERE/'candidates/paired_chemistry'


def digest(p):
    with p.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    build = json.loads((HERE/'results/paired_chemistry_build.json').read_text())
    assert build['physical_npz_unchanged']
    assert all(digest(Path(p))==sha for p,sha in build['source_sha256'].items())
    files = {p.name:digest(p) for p in DIRECTORY.iterdir() if p.is_file()}
    words = np.load(DIRECTORY/'submission.npz',allow_pickle=False)['angle_words']
    circuits = load_circuits(DIRECTORY)
    bit = lambda a,n: ((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    rows = []
    for seed in [2026102533,2026102534]:
        rng = np.random.default_rng(seed)
        address = np.tile(np.arange(len(words)),2)
        for k,(circuit,delta) in enumerate(zip(circuits,batch_deltas(words,4))):
            assert digest(DIRECTORY/f'angle_k_{k}.kmb')==build['rows'][k]['sha256']
            sim = km.Simulator(batch_size=len(address),seed=seed+k)
            sim.use_same_registers_as(circuit)
            values,public = {},set()
            for name,reg in circuit.register_data.items():
                public.update(reg)
                values[name] = address if name=='selection' else rng.integers(0,256,len(address))
                sim.write_across_shots(reg,bit(values[name],len(reg)))
            scratch = km.array([km.q(i) for i in range(circuit.num_qubits) if km.q(i) not in public])
            for reverse in [False,True]:
                sim.do(circuit)
                for name,reg in circuit.register_data.items():
                    want = values[name] if name=='selection' or reverse else values[name]^delta[address,int(name[5:])]
                    assert np.array_equal(sim.read_across_shots(reg),bit(want,len(reg))),(seed,k,name,reverse)
                assert not sim.read_across_shots(scratch).any(),(seed,k,'scratch',reverse)
                assert not sim.read_phase_flipped_across_shots().any(),(seed,k,'phase',reverse)
            rows.append(dict(seed=seed,k=k,shots_per_direction=len(address),status='PASS'))
    (HERE/'results/paired_chemistry_components.json').write_text(json.dumps(dict(status='PASS',rows=rows,all_valid_addresses=True,arbitrary_targets=True,paired_application=True,phase_and_scratch=True),indent=2)+'\n')
    print(json.dumps(dict(stage='ALL_ADDRESS_DIRTY_TARGET_COMPONENT_CHECKS_PASS',circuits=len(circuits),seeds=2)),flush=True)
    verifier = (HERE/'../../../from_host/verifiers/dfthcsos/verify.py').resolve()
    host_sha = digest(verifier)
    runs = []
    for run,salt in [('default',''),('independent','agent20-paired-chemistry-independent-20261025')]:
        command = [sys.executable,str(verifier),str(DIRECTORY)]
        if salt:
            command += ['--salt',salt]
        print(json.dumps(dict(stage='FULL_NATIVE_START',run=run,salt=salt,files=files,verifier_sha256=host_sha)),flush=True)
        output = HERE/'results'/f'paired_chemistry_{run}_verify.json'
        start = time.time()
        with output.open('w') as out,output.with_suffix('.err').open('w') as err:
            result = subprocess.run(command,stdout=out,stderr=err,cwd=HERE)
        report = json.loads(output.read_text())
        row = dict(run=run,salt=salt,seconds=time.time()-start,returncode=result.returncode,
                   status=report['status'],score=report.get('score'),submission_sha256=report.get('submission_sha256'))
        runs.append(row)
        (HERE/'results/paired_chemistry_verification_runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        print(json.dumps(row),flush=True)
        assert result.returncode==0 and report['status']=='PASS' and report['verified_shots']==50880
        assert report['toffoli_count']==build['total_Toffolis'] and report['num_qubits']==454
        assert host_sha==digest(verifier) and all(digest(DIRECTORY/p)==sha for p,sha in files.items())
    assert runs[0]['submission_sha256']==runs[1]['submission_sha256']
    print(json.dumps(dict(status='TWO_FULL_NATIVE_CHECKS_COMPLETE',runs=runs)),flush=True)


if __name__=='__main__':
    main()
