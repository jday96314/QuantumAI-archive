"""Rebuild the optimized submission from the frozen parameters and gate builders."""
import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'from_host/kickmix'),str(ROOT/'from_host/verifiers'),str(HERE)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=HERE/'submission')
    p.add_argument('--task',choices=['all','invsqrt','dqi','dfthcsos'],default='all')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    params=json.loads((HERE/'parameters/solution.json').read_text())
    import kickmix as km
    import numpy as np
    if a.task in ['all','invsqrt']:
        import polynomial_invsqrt as poly
        settings=params['invsqrt'];widths=tuple(settings['widths']);frac=settings['frac']
        table=np.load(HERE/'parameters/invsqrt_table.npy',allow_pickle=False)
        original=poly.table_for
        def fitted(f,w,b=0):
            if f==frac and tuple(w)==widths:return table
            return original(f,w,b)
        poly.table_for=fitted
        try:circuit=poly.build(frac,widths,0)
        finally:poly.table_for=original
        circuit.to_file(a.output/'invsqrt.kmb',format='kmb')
        print('invsqrt',circuit.num_qubits,circuit.operation_statistics().num_toffoli_like,flush=True)
    if a.task in ['all','dqi']:
        from fast_dqi import FastBuilder
        if params['dqi'].get('forney'):
            from forney_dqi import ForneyDecoder as Decoder
        elif params['dqi']['compact']:
            from compact_dqi import CompactDecoder as Decoder
        else:
            from dqi.dqi_reference_circuit import ProosAndZalkaDecoder as Decoder
        b=FastBuilder(km.CircuitBuilder())
        syndromes=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
        mbuc=[b.create_classical_register(12,f'E{i}') for i in range(4095)]
        Decoder(70,km.GF2Field(12)).build_syndrome_decoding_mbuc(b,syndromes=syndromes,error_mbuc=mbuc,n_iters=params['dqi']['iterations'])
        circuit=b.finish_circuit();circuit.to_file(a.output/'dqi.kmb',format='kmb')
        print('dqi',circuit.num_qubits,circuit.operation_statistics().num_toffoli_like,flush=True)
    if a.task in ['all','dfthcsos']:
        from chemistry_qrom import build
        from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
        path=a.output/'dfthcsos';path.mkdir(exist_ok=True)
        # Remove only prior generated angle-family files, avoiding a mixed family.
        for old in path.glob('angle_k_*.kmb'):old.unlink()
        shutil.copyfile(HERE/'parameters/chemistry.npz',path/'submission.npz')
        with np.load(path/'submission.npz',allow_pickle=False) as data:words=data['angle_words']
        for i,delta in enumerate(batch_deltas(words,params['dfthcsos']['batch_size'])):
            circuit=build(delta,params['dfthcsos']['brot'])
            assert circuit.num_qubits<=128
            circuit.to_file(path/f'angle_k_{i}.kmb',format='kmb')
            print('dfthcsos',i,circuit.num_qubits,circuit.operation_statistics().num_toffoli_like,flush=True)

if __name__=='__main__':main()
