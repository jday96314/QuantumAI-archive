"""Rebuild frozen research candidates, or run unchanged full host verification."""
import argparse,hashlib,json,subprocess,sys
from model import HERE
from build import build
from optimize_workspace import compact

def main():
    p=argparse.ArgumentParser();p.add_argument('--full',action='store_true');a=p.parse_args();checks={}
    for mode,name in (('split','split_q84_k0'),('rows','rows_q84_k0')):
        frozen=HERE/'candidates'/f'{name}.kmb';rebuilt=HERE/'candidates'/f'{name}_rebuilt.kmb'
        c,_=build(mode,84,0);c=compact(c,rebuilt);c.to_file(str(rebuilt),format='kmb')
        assert frozen.read_bytes()==rebuilt.read_bytes(),name
        checks[name]=hashlib.sha256(frozen.read_bytes()).hexdigest();rebuilt.unlink()
        if a.full:
            out=HERE/'results'/f'{mode}_full_verify.json'
            with out.open('w') as f:
                subprocess.run([sys.executable,str(HERE.parents[2]/'from_host/verifiers/invsqrt/verify_submission_invsqrt.py'),
                    str(frozen),'--max-workers','2','--salt',f'inv-phase-{mode}-final-20261024'],stdout=f,check=True)
    (HERE/'results/rebuild.json').write_text(json.dumps(dict(status='PASS',byte_identical=checks),indent=2)+'\n')
    print(json.dumps(checks))

if __name__=='__main__':main()
