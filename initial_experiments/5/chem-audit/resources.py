"""Host-exact PREPARE/Givens and constructive ANGLE upper bounds.

Reuses chem-blocks/signed_oracle.py, read-only. No physical feasibility is
inferred from this resource certificate.
"""
import importlib.util,json,math
from audit_core import *
from dfthcsos.verify import heuristic_costs
spec=importlib.util.spec_from_file_location('audit_signed_bound',HERE.parent/'chem-blocks/signed_oracle.py')
oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
rows=[]
for R,B,C in [(1,719,350),(1,639,350),(2,447,175),(2,463,175),(2,479,175),(4,180,88),(1,511,350)]:
    for bits in [7,8]:
        h=heuristic_costs(76,R,B,C,bits);options=[]
        for J in range(4,12):
            ob=oracle.bound(76+R*B,bits,J)
            if ob is None:continue
            T=h['toffoli']+2*math.ceil(75/J)*ob['T_per_oracle']
            Q=h['persistent_qubits']+J*bits+max(h['temporary_qubits'],ob['angle_temp'])
            options.append(dict(R=R,B=B,C=C,bits=bits,J=J,T=T,Q=Q,oracle=ob,host=h,
                conditional_cost_incumbent_lambda=T*Q*48.59417497035283,
                conditional_cost_at_normalization_bound=T*Q*46.78617959252637,
                scope='host-exact PREPARE/Givens plus realizable oracle upper bound; physical feasibility untested'))
        rows.append(min(options,key=lambda x:x['T']*x['Q']))
write_json(HERE/'results/resources.json',rows)
for r in rows:print(json.dumps({k:v for k,v in r.items() if k not in ['host','oracle']}),flush=True)
