"""Freeze verified complete submissions at thresholds 10 through 14."""
import hashlib,json,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
 summary=json.loads((HERE/'submission.summary.json').read_text());audit=json.loads((HERE/'results/archive_audit.json').read_text())
 assert summary['status']==audit['status']=='PASS'
 assert summary['sha256']==audit['archive_sha256']==sha(HERE/'submission.zip')
 assert summary['overall_score']==audit['overall_score']
 assert summary['verified_shots']==dict(invsqrt=200540928,dqi=102400,dfthcsos=50880)
 saved=[]
 for threshold in range(10,15):
  dest=HERE/'milestones'/str(threshold)
  if summary['overall_score']<threshold or dest.exists():continue
  dest.mkdir(parents=True)
  shutil.copyfile(HERE/'submission.zip',dest/'submission.zip')
  assert sha(dest/'submission.zip')==summary['sha256']
  for name,source in [('submission.summary.json',HERE/'submission.summary.json'),('archive_audit.json',HERE/'results/archive_audit.json')]:shutil.copyfile(source,dest/name)
  (dest/'milestone.json').write_text(json.dumps(dict(threshold=threshold,score=summary['overall_score'],archive_sha256=summary['sha256']),indent=2)+'\n')
  saved.append(threshold)
 result=dict(best_verified_local_score=summary['overall_score'],new_milestones=saved,remaining_thresholds=[i for i in range(10,15) if not (HERE/'milestones'/str(i)/'submission.zip').exists()])
 (HERE/'milestones/progress.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
