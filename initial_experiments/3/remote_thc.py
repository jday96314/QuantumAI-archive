"""Read selected files from Lee et al.'s CC-BY-4.0 THC data via HTTP ranges.

Source: https://doi.org/10.5281/zenodo.4248322
Only data are downloaded; no code from the archive is executed.
"""
import argparse,io,json,tarfile
from pathlib import Path
import requests
HERE=Path(__file__).resolve().parent
URL='https://zenodo.org/api/records/4248322/files/qc_thc.tar/content'
SIZE=9211981312
class Remote(io.RawIOBase):
 def __init__(self):self.pos=0;self.offset=-1;self.data=b'';self.session=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else SIZE+offset
  return self.pos
 def read(self,n=-1):
  if n<0:n=SIZE-self.pos
  if n>100_000_000:raise ValueError('explicit selection required for large downloads')
  result=[]
  while n and self.pos<SIZE:
   if not(self.offset<=self.pos<self.offset+len(self.data)):
    start=self.pos;end=min(SIZE-1,start+max(16384,min(n,8_000_000))-1)
    r=self.session.get(URL,headers={'Range':f'bytes={start}-{end}'},timeout=60)
    r.raise_for_status()
    if r.status_code!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {start}-'):
     raise ValueError(('server did not honor range',r.status_code,r.headers.get('Content-Range')))
    self.offset=start;self.data=r.content
   k=min(n,self.offset+len(self.data)-self.pos)
   result.append(self.data[self.pos-self.offset:self.pos-self.offset+k]);self.pos+=k;n-=k
  return b''.join(result)
def main():
 p=argparse.ArgumentParser();p.add_argument('--extract',nargs='*');a=p.parse_args();root=HERE/'vendor/published_thc'
 root.mkdir(parents=True,exist_ok=True);index=[]
 with tarfile.open(fileobj=Remote(),mode='r:') as archive:
  for entry in archive:
   if entry.isfile():
    index.append(dict(name=entry.name,size=entry.size,offset=entry.offset_data))
    print(json.dumps(index[-1]),flush=True)
    if a.extract and entry.name in a.extract:
     name=Path(entry.name).name
     if name.startswith('._'):continue
     data=archive.extractfile(entry).read();(root/name).write_bytes(data)
 (root/'index.json').write_text(json.dumps(index,indent=2)+'\n')
if __name__=='__main__':main()
