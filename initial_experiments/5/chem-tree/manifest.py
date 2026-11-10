"""Hash delivered project artifacts; caches/temp and this manifest are excluded."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
excluded={'tmp','.cache','.cuda-cache','.matplotlib','__pycache__'}
rows=[]
for p in sorted(HERE.rglob('*')):
    relative=p.relative_to(HERE)
    if not p.is_file() or p.name=='artifact_manifest.json' or any(x in excluded for x in relative.parts):continue
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    rows.append(dict(path=str(relative),bytes=p.stat().st_size,sha256=h.hexdigest(),
                quarantined='aborted_encoder' in relative.parts))
(HERE/'artifact_manifest.json').write_text(json.dumps(dict(artifacts=rows,exclusions=sorted(excluded),
        note='No self hash; immutable external sources are bound in results/reference_provenance.json'),indent=2)+'\n')
print(json.dumps(dict(artifacts=len(rows),bytes=sum(r['bytes'] for r in rows))))
