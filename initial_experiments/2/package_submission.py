"""Package the fully verified submission, checking archive CRCs and hashes."""
import hashlib
import json
from pathlib import Path
import zipfile
HERE=Path(__file__).resolve().parent

def main():
    validation=json.loads((HERE/'results/final-validation/summary.json').read_text())
    contract=json.loads((HERE/'results/metric-contract.json').read_text())
    assert validation['status']==contract['status']=='PASS'
    archive=HERE/'submission.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for path in sorted((HERE/'submission').rglob('*')):
            if path.is_file():z.write(path,path.relative_to(HERE/'submission').as_posix())
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert set(z.namelist())==set(contract['files'])
    digest=hashlib.sha256()
    with archive.open('rb') as f:
        for chunk in iter(lambda:f.read(16*1024*1024),b''):digest.update(chunk)
    result={'archive':'submission.zip','bytes':archive.stat().st_size,'sha256':digest.hexdigest(),'local_overall_score':validation['overall_score'],'status':'PASS'}
    (HERE/'results/archive.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
