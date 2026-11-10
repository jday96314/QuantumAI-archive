"""Compile complete actual ANGLE tables and host costs for physical screens.

Invalid physical models remain invalid: these are exact static resource
products, not accepted chemistry scores or new native circuit PASS claims.
"""
from dataclasses import replace
from core import *
from signed_oracle import build,bound
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.verify import resource_costs

def main():
    inp=load_inputs(INC);records=[]
    for path in sorted((HERE/'parameters').glob('screen_*_physical.npz')):
        tag=path.stem.removesuffix('_physical');z=dict(np.load(path));physical=json.loads((HERE/f'results/{tag}_physical.json').read_text())
        params={k:z[k] for k in ['U','W','H2_I','H1_I','Hsym_lower']};w=z['W'];off=z['H2_I']-w.sum(1)
        this=replace(inp,params=params,w=w,W=off);A=len(z['angle_words']);ob=bound(A,8,7);circuits=[]
        out=HERE/f'resource_circuits/{tag}';out.mkdir(parents=True,exist_ok=True)
        for k,delta in enumerate(batch_deltas(z['angle_words'],7)):
            c=build(delta,8,ob['lowbits']);assert c.num_qubits<=128
            c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');circuits.append(c)
        T=2*sum(int(c.operation_statistics().num_toffoli_like) for c in circuits)
        resources=resource_costs(this,circuits,8,T)
        rec=dict(physical_screen=tag,physical_pass=physical['passes_physical_gates'],
                 resources_static=resources,lambda_eff=physical['native_metrics']['lambda_eff'],
                 invalid_resource_product=physical['native_metrics']['lambda_eff']*resources['toffoli']*resources['qubits'],
                 accepted_task_score=None,circuits=str(out),native_family_verification=False)
        records.append(rec);print(json.dumps(rec),flush=True)
    report(HERE/'results/actual_resources.json',dict(scope='Full actual ANGLE families compiled; exact static host resources; physical failures remain unscored',records=records))
if __name__=='__main__':main()
