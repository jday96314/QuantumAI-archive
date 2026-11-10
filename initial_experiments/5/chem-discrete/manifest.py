"""Final artifact hashes; only this experiment directory is traversed."""
import hashlib,json
from pathlib import Path
here=Path(__file__).resolve().parent
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
records=[]
excluded=['results/plot_dependencies/','results/mplconfig/','results/pip_tmp/']
scope=dict(description='Research artifacts; regenerable plotting dependencies/cache and the manifest stdout log are excluded.',excluded_prefixes=excluded,excluded_files=['artifact_manifest.json','results/manifest.log'])
(here/'manifest_scope.json').write_text(json.dumps(scope,indent=2))
for p in sorted(here.rglob('*')):
    rel=str(p.relative_to(here))
    if p.is_file() and rel not in scope['excluded_files'] and not any(rel.startswith(x) for x in excluded):records.append(dict(path=rel,bytes=p.stat().st_size,sha256=digest(p)))
(here/'artifact_manifest.json').write_text(json.dumps(records,indent=2));print(json.dumps(dict(files=len(records),bytes=sum(r['bytes'] for r in records))))
