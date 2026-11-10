"""Exact emitted gate/lifetime profiles without retaining the operation tape."""
import argparse,json
from pathlib import Path
import hybrid_decoder as decoder
from import_published import load


class DiscardTape:
    def __init__(self):self.n=0
    def __len__(self):return self.n
    def _append_row(self,*args):self.n+=1
    def append(self,*args):self.n+=1


def module():
    m=load();Builder=m.Builder
    class CountBuilder(Builder):
        def __init__(self):super().__init__();self._items=DiscardTape()
    m.Builder=CountBuilder;return m


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--after',type=int,default=1);p.add_argument('--block',type=int,default=3);p.add_argument('--tag',default='');a=p.parse_args()
    decoder.load=module;decoder.build.profile_only=True
    r=decoder.build(70,a.block,a.after,4095,Path(f'candidates/profile_b{a.block}_after{a.after}{a.tag}.kmb'),True,True,True,True,True)
    Path(f'results/profile_b{a.block}_after{a.after}{a.tag}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='rows'}),flush=True)
