"""Resource breakdown and evidence binding; writes only experiment-local files."""
import csv,hashlib,json,math
from pathlib import Path
import kickmix as km
from model import HERE,synthesis

def stages(path):
    c=km.Circuit.from_file(str(path));tmp=HERE/'results'/f'{path.stem}_profile.kmx';c.to_file(str(tmp),format='kmx')
    marks=set(json.loads((HERE/'results'/f'{path.stem}_build.json').read_text())['profile']['peaks'])
    out={};current='start'
    for line in tmp.read_text().splitlines():
        if line.startswith('INFO_TAG '):
            mark=json.loads(line[len('INFO_TAG '):])
            if mark in marks:current=mark
        op=line.split()[0] if line else ''
        if op in ('CCX','CCZ','CSWAP'):out[current]=out.get(current,0)+1
    assert sum(out.values())==c.operation_statistics().num_toffoli_like
    tmp.unlink();return out

def main():
    paths=['split_q84_k0','rows_q84_k0'];rows=[];audit={}
    for name in paths:
        path=HERE/'candidates'/f'{name}.kmb';record=json.loads((HERE/'results'/f'{name}_build.json').read_text())
        record['static_Toffolis_by_stage']=stages(path)
        (HERE/'results'/f'{name}_build.json').write_text(json.dumps(record,indent=2)+'\n')
        full=HERE/'results'/('split_full_verify.json' if name.startswith('split') else 'rows_full_verify.json')
        if full.exists() and full.stat().st_size:
            result=json.loads(full.read_text());assert result['status']=='PASS'
            assert result['verified_shots']==200540928
            assert result['circuit_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
            model=json.loads((HERE/'results/model.json').read_text())
            assert result['worst_error_over_tolerance']==model['worst_error_over_tolerance']
            assert abs(result['rotation_eps']*result['rotation_count']*2**64/math.pi-model['min_margin_units'])<1e-6
            assert abs(result['score']-result['num_qubits']*result['toffoli_equivalent'])<1e-6
            row=dict(name=name,Q=result['num_qubits'],static_Toffolis=result['max_toffoli_count'],
                mean_Toffolis=result['toffoli_count'],rotations=result['rotation_count'],
                synthesis_work=(result['toffoli_equivalent']-result['toffoli_count']),cost=result['score'],native_shots=result['verified_shots'])
        else:
            row=dict(name=name,Q=record['Q'],static_Toffolis=record['static_Toffolis'],
                mean_Toffolis=None,rotations=record['rotation_count'],synthesis_work=record['synthesis_work'],cost=record['estimated_cost'],native_shots=0)
        rows.append(row);audit[name]=dict(candidate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),native_result=str(full),resources=row)
    fallback=HERE/'reference/incumbent.kmb';sha=hashlib.sha256(fallback.read_bytes()).hexdigest()
    assert sha=='27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1'
    fallback_result=json.loads((HERE/'reference/incumbent_verify.json').read_text())
    assert fallback_result['circuit_sha256']==sha
    assert fallback_result['verified_shots']==200540928 and fallback_result['status']=='PASS'
    audit['fallback_sha256']=sha
    dependencies=json.loads((HERE/'results/read_only_dependencies.json').read_text())
    assert all(hashlib.sha256((HERE.parents[2]/p).read_bytes()).hexdigest()==h for p,h in dependencies.items())
    audit['read_only_dependencies']=dependencies
    audit['rebuild']=json.loads((HERE/'results/rebuild.json').read_text())
    assert all(audit['rebuild']['byte_identical'][name]==audit[name]['candidate_sha256'] for name in paths)
    audit['total_experimental_full_shots']=sum(r['native_shots'] for r in rows)
    audit['status']='PASS' if audit['total_experimental_full_shots']==401081856 else 'INCOMPLETE'
    audit['sources']={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')}
    audit['references']={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'reference').glob('*') if p.suffix!='.npz'}
    (HERE/'results/audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    with (HERE/'results/resource_frontier.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(resources=rows,stages={name:stages(HERE/'candidates'/f'{name}.kmb') for name in paths})))

if __name__=='__main__':main()
