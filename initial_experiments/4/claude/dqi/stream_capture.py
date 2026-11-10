"""Syndrome-consuming adaptation of the attributed upstream q-form BM.

The register contract stays unchanged. Previously consumed syndromes become
backup-bank storage. A single vanished leading coefficient parks one input;
later exceptional paths stop with a reversible flag.
"""
import json,time
from contextlib import contextmanager
import numpy as np
from import_published import load,HERE

def capture(output='stream'):
    m=load();OldBank=m.BankReg;OldBM=m.BM
    compact='compact' in output
    bank_offset=3 if 'offset3' in output else 4 if 'offset4' in output else 5
    consume_limit={3:67,4:65,5:63}[bank_offset]
    old_swap_fn=m.Em.swap_fn
    def frozen_swap_fn(self,control,pairs):
        return old_swap_fn(self,control,tuple(pairs))
    m.Em.swap_fn=frozen_swap_fn
    class Alloc(m.Alloc):
        def __init__(self,b,syn):
            super().__init__(b);self.syn=syn;self.bb_count=0;self.borrowed=set()
        def close(self,name,qs):
            if set(qs).issubset(self.borrowed):return
            super().close(name,qs)
    class Bank(OldBank):
        __slots__=('slot',)
        def __init__(self,A,name,width=m.B):
            super().__init__(A,name,width);self.slot=None
            if name=='bb':self.slot=A.bb_count;A.bb_count+=1
        def mat(self):
            if self.wires is None and self.slot is not None and self.slot>=bank_offset:
                A=self.A;self.wires=list(A.syn[self.slot-bank_offset]);A.borrowed.update(self.wires)
                self.round=A.round;A.bank_log.setdefault(A.round,[]).append(self)
                return self.wires
            return super().mat()
    m.BankReg=Bank
    class BM(OldBM):
        def __init__(self,E,A,syn,Dw,dW):
            super().__init__(E,A,syn,Dw,dW)
            self.active=E.b.qubit();self.filled=E.b.qubit()
            self.stopidx=E.b.qubits(9);self.parkidx=E.b.qubits(6);self.park=E.b.qubits(12)
            self.held=[self.active,self.filled,*self.stopidx,*self.parkidx,*self.park]
        @contextmanager
        def tmp(self,name,n):
            with m.opened(self.E,name,n) as qs:yield qs
        @contextmanager
        def inverse_live(self,source,target):
            E=self.E
            if not compact:
                with m.chain_live(E,source,target):yield
                return
            wins,rts=m.chain_fwd(E,source,target)
            E.frees(rts[:-1]);m.chain_close(E,wins)
            yield
            for (name,_),window in zip(m.CHAIN_WINS,wins):E.A.reopen(name,window)
            E.macs(rts[:-1]);E.free(rts[-1]);E.frees(rts[:-1]);m.chain_close(E,wins)
        def rt(self,a,b,target,ca=0,cb=0):
            return [m.norm_term(tuple(a),ca,tuple(b),cb,tuple(target))]
        def stop_bad(self,bad,code):
            E=self.E
            E.cx(bad,self.active)
            for j in m.bitlist(code,9):E.cx(bad,self.stopidx[j])
            with m.eq_test(E,self.stopidx,code) as flag:E.free([('L',(flag,),0,(bad,))])
        def bad_nonzero(self,word,bad,gate):
            with m.eq_test(self.E,list(word),0) as zero:
                self.E.mac(self.rt((zero,),(gate,),(bad,),ca=1))
        def stop_word(self,word,code,sign=None):
            E=self.E
            with self.tmp('sbad',1) as (bad,):
                if sign is None:self.bad_nonzero(word,bad,self.active)
                else:
                    with self.tmp('sgate',1) as (gate,):
                        rt=self.rt((self.active,),(sign,),(gate,));E.mac(rt)
                        self.bad_nonzero(word,bad,gate);E.free(rt)
                self.stop_bad(bad,code)
        def round_head(self,r,bb):
            A,E=self.A,self.E;A.round=r;A.log=self.pools.setdefault(r,set())
            self.JROUND.append(self.J.mark());self.STK[0].mat()
            size=min(37,max(m.KMSV[r],m.KMMV[r])+2)
            deferred=bank_offset==3 and r%2 and r<=consume_limit
            for k in range(size):
                if not (deferred and bb[k].wires is None and bb[k].slot>=bank_offset):bb[k].mat()
            if compact:
                for k in range(max(m.IMAXV[r],m.KMSV[r],m.KMMV[r])+1):self.SIG[k].mat()
            if deferred:
                if r and r<69:
                    if bb[size-1].wires is not None:self.stop_word(bb[size-1],4*r+1)
                    self.stop_word(self.STK[0],4*r+2,sign=self.dW[m.SGN])
                return bb
            self.rotate_bank(r,bb)
            if r and r<69:
                self.stop_word(bb[0],4*r+1)
                self.stop_word(self.STK[0],4*r+2,sign=self.dW[m.SGN])
            return bb
        def rotate_bank(self,r,bb):
            size=min(37,max(m.KMSV[r],m.KMMV[r])+2)
            for k in range(size):bb[k].mat()
            for k in reversed(range(1,size)):self.E.swaps(None,zip(bb[k],bb[k-1]))
        def discrepancy(self,r):
            E=self.E
            def plain():
                if 'gated' not in output:
                    super(BM,self).discrepancy(r)
                    return
                # After completion, the tighter bank can reach an unconsumed
                # input. Freeze the discrepancy as well as the update on those
                # branches, so the bank never spoils a later erase claim.
                with self.tmp('rawD',12) as raw:
                    rt=[]
                    for i in range(m.IMAXV[r]+1):rt+=m.prod_rt(raw,self.SIG[i],self.syn[r-i])
                    E.mac(rt)
                    E.mac([m.norm_term((self.active,),0,(raw[j],),0,(self.Dw[j],)) for j in range(12)])
                    E.free(rt)
            if r%2==0 or r>consume_limit:
                # Only the even round immediately following a parked input
                # still needs that original syndrome.
                if r%2==0 and r>0:
                    with m.eq_test(E,self.parkidx,r//2) as equal:
                        with self.tmp('restoregate',1) as (gate,):
                            rt=self.rt((equal,),(self.active,),(gate,));E.mac(rt)
                            E.swaps(gate,zip(self.syn[r//2-1],self.park))
                            plain()
                            E.swaps(gate,zip(self.syn[r//2-1],self.park));E.free(rt)
                else:plain()
                return
            l=(r+1)//2;target=self.syn[l-1]
            # A second parking request stops before changing an input.
            with self.tmp('sbad',1) as (bad,):
                with self.tmp('sgate',1) as (gate,):
                    rt=self.rt((self.active,),(self.filled,),(gate,));E.mac(rt)
                    with m.eq_test(E,list(self.SIG[l]),0) as zero:
                        E.mac(self.rt((zero,),(gate,),(bad,)))
                    E.free(rt)
                self.stop_bad(bad,4*r+3)
            with self.tmp('parkctl',1) as (parkctl,):
                with m.eq_test(E,list(self.SIG[l]),0) as zero:
                    E.mac(self.rt((zero,),(self.active,),(parkctl,)))
                E.swaps(parkctl,zip(target,self.park));E.cx(parkctl,self.filled)
                for j in m.bitlist(l,6):E.cx(parkctl,self.parkidx[j])
                with m.eq_test(E,self.parkidx,l) as equal:
                    E.free([('L',(equal,),0,(parkctl,))])
            # Arm the variable diagonal against zero and inactive branches.
            with m.eq_test(E,list(self.SIG[l]),0) as zero:
                E.cx(zero,self.SIG[l][0])
                if 'noscale' in output:
                    scale=self.SIG[l]
                    with self.tmp('ssinv',12) as iv,m.chain_live(E,scale,iv),self.tmp('ssproduct',12) as product:
                        def forward():m.mon_slot(E,target,scale,iv,product)
                        def backward():m.mon_slot(E,target,iv,scale,product)
                        E.J.step(forward,backward)
                else:
                    with self.tmp('sscale',12) as scale:
                        rt=[m.norm_term((self.active,),0,(self.SIG[l][j],),int(j==0),(scale[j],)) for j in range(12)]
                        E.mac(rt);E.x(scale[0])
                        with self.tmp('ssinv',12) as iv,m.chain_live(E,scale,iv),self.tmp('ssproduct',12) as product:
                            def forward():m.mon_slot(E,target,scale,iv,product)
                            def backward():m.mon_slot(E,target,iv,scale,product)
                            E.J.step(forward,backward)
                        E.x(scale[0]);E.free(rt)
                E.cx(zero,self.SIG[l][0])
            with self.tmp('sgatedsigma',12) as gated:
                for i in range(m.IMAXV[r]+1):
                    if i==l:continue
                    rt=[m.norm_term((self.active,),0,(self.SIG[i][j],),0,(gated[j],)) for j in range(12)]
                    E.mac(rt);E.xor(m.prod_rt(target,gated,self.syn[r-i]));E.free(rt)
            E.swaps(self.active,zip(target,self.Dw))
            # The first zero diagonal with a zero discrepancy is completion.
            with self.tmp('sbad',1) as (bad,):
                with m.eq_test(E,self.parkidx,l) as equal:
                    with self.tmp('sgate',1) as (gate,):
                        rt=self.rt((equal,),(self.active,),(gate,));E.mac(rt)
                        with m.eq_test(E,self.Dw,0) as zero:
                            E.mac(self.rt((zero,),(gate,),(bad,)))
                        E.free(rt)
                self.stop_bad(bad,4*r+4)
        @contextmanager
        def quotient(self,ctl,bb):
            E=self.E
            with self.tmp('pq',12) as pv:
                E.mac(m.prod_rt(pv,self.Dw,self.SIG[0]))
                if 'select' in output:
                    # Select at the six-bit inverse norm, before the final
                    # field product. ctl implies active, so the other branch
                    # is simply (active XOR ctl) * pv.
                    with self.tmp('m',12) as mq:
                        wins,rts=m.chain_open(E,pv,mq);E.macs(rts[:-1])
                        with self.tmp('selectednorm',6) as gy:
                            rg=[m.norm_term((ctl,),0,(wins[-1][j],),0,(gy[j],)) for j in range(6)]
                            E.mac(rg)
                            rt=m.chain_rts(pv,[*wins[:-1],gy],mq)[-1]
                            rt += [m.norm_term(tuple(sorted((self.active,ctl))),0,(pv[j],),0,(mq[j],)) for j in range(12)]
                            E.mac(rt);yield mq;E.free(rt);E.free(rg)
                        E.frees(rts[:-1]);m.chain_close(E,wins)
                else:
                    with self.tmp('dinv',12) as dv,self.inverse_live(pv,dv),self.tmp('mraw',12) as raw,self.tmp('m',12) as mq:
                        rt=[]
                        for j in range(12):
                            rt.append(('L',(pv[j],),0,(raw[j],)))
                            rt.append(m.norm_term((ctl,),0,tuple(sorted((pv[j],dv[j]))),0,(raw[j],)))
                        E.mac(rt)
                        gated=[m.norm_term((self.active,),0,(raw[j],),0,(mq[j],)) for j in range(12)]
                        E.mac(gated);yield mq;E.free(gated);E.free(rt)
                with self.tmp('gatedbb',12) as gb:
                    rtg=[m.norm_term((self.active,),0,(bb[0][j],),0,(gb[j],)) for j in range(12)]
                    E.mac(rtg)
                    rt=m.prod_rt(pv,self.Dw,self.SIG[0])+m.prod_rt(pv,self.Dw,gb)
                    rt.append(('L',(ctl,),0,(pv[0],)))
                    E.free(rt);E.free(rtg)
        def push(self,r,ctl=True):
            # The pre-round guard guarantees an empty slot on active pushes.
            E=self.E
            with self.tmp('pushgate',1) as (gate,):
                rt=(self.rt((self.active,),(self.dW[m.SGN],),(gate,)) if ctl else
                    [('L',(self.active,),0,(gate,))])
                E.mac(rt)
                if m.CAPS[r]:
                    E.swaps(gate,zip(self.Dw,self.STK[0][:12]))
                    if r<69:
                        for j in m.bitlist(r,7):E.cx(gate,self.STK[0][12+j])
                E.free(rt)
        def close_ctl(self,r,ctl):
            E=self.E
            for i in range(m.DWB):E.cx(ctl,self.dW[i])
            if m.CAPS_FULL[r]==0:
                E.free(self.rt((self.active,),(self.dW[m.SGN],),(ctl,)))
            else:
                with self.tmp('closegate',1) as (gate,):
                    rt=self.rt((self.active,),(self.dW[m.SGN],),(gate,));E.mac(rt)
                    with m.eq_test(E,list(self.STK[0][12:]),r) as equal:
                        E.free(self.rt((gate,),(equal,),(ctl,),cb=1))
                    E.free(rt)
        def round(self,r,bb):
            E=self.E
            with self.b.phase(f'bm/r{r:02d}'):
                bb=self.round_head(r,bb)
                with self.tmp('ctl',1) as (ctl,):
                    self.discrepancy(r)
                    if bank_offset==3 and r%2 and r<=consume_limit:self.rotate_bank(r,bb)
                    with self.tmp('nonzero',1) as (nz,),self.tmp('ctlguard',1) as (guard,):
                        m.flag_test(E,self.Dw,(),nz)
                        rt=self.rt((self.active,),(self.dW[m.SGN],),(guard,),cb=1);E.mac(rt)
                        E.mac(self.rt((nz,),(guard,),(ctl,)));E.free(rt)
                        # Bare nonzero flag is a pure function of Dw.
                        with m.eq_test(E,self.Dw,0) as zero:E.free([('L',(zero,),1,(nz,))])
                    with self.quotient(ctl,bb) as mq:
                        self.bank_swaps(r,ctl,bb);self.update(r,0,mq,bb)
                    self.push(r);self.erase_dw(ctl);self.close_ctl(r,ctl)
                if r<68:self.inc_d()
            return bb
        def round_last(self,bb):
            E=self.E;r=69
            with self.b.phase('bm/r69'):
                bb=self.round_head(r,bb);self.discrepancy(r)
                with self.tmp('lastp',12) as pv,self.tmp('m',12) as mq:
                    rt=m.prod_rt(pv,self.Dw,self.SIG[0]);E.mac(rt)
                    rg=[m.norm_term((self.active,),0,(pv[j],),0,(mq[j],)) for j in range(12)]
                    E.mac(rg);self.update(r,1,mq,bb);E.free(rg);E.free(rt)
                self.push(r,ctl=False)
        def run(self):
            E,J,A=self.E,self.J,self.A;J.recording=True
            A.round=-1;A.log=self.pools.setdefault(-1,set())
            with E.using(J=J,A=A):
                with self.b.phase('bm/seed'):
                    E.x(self.SIG[0][0]);E.x(self.BB[0][0]);E.x(self.active)
                bb=self.round0(list(self.BB))
                for r in range(1,69):bb=self.round(r,bb)
                self.round_last(bb)
            J.recording=False;self.JROUND.append(J.mark());A.round=None;A.log=None
    start=time.time();b=m.Builder();marks=[];old_phase=b.phase;PEAKS=[];TOFS=[]
    @contextmanager
    def marked_phase(name):
        marks.append((len(b._items),name))
        TOFS.append((name,b._stream_totals()[0]))
        before=b.active_qubits
        with old_phase(name):yield
        PEAKS.append((name,before,b.active_qubits))
    b.phase=marked_phase
    syn=list(m.declare_registers(b,[(f's{i:02d}',12) for i in range(70)]).values())
    A=Alloc(b,syn);held=b.qubits(23);Dw=held[:12];dW=held[12:20]
    E=m.Em(b,*held[20:23],A)
    bm=BM(E,A,syn,Dw,dW);bm.run()

    import collections
    owner={};live=set();snap=None;kinds=collections.Counter();allocat={};curmark=None
    for ev in b.wire_log:
        kind=ev[0];kinds[kind]+=1
        if kind=='mark':curmark=ev[1]
        if kind=='alloc':
            live.add(ev[1])
            allocat[ev[1]]=curmark
            if ev[3]==b.peak_qubits and snap is None:
                snap=collections.Counter(owner.get(q,'?').split('/')[0] for q in live)
                synq=set(q for w in syn for q in w)|set(held)
                unk=[q for q in live if q not in owner and q not in synq]
                print('UNOWNED',len(unk),collections.Counter(allocat.get(q) for q in unk),flush=True)
        elif kind in('free','recycle'):
            live.discard(ev[1]);owner.pop(ev[1],None)
        elif kind=='claim':
            for q in ev[1][1]:owner[q]=ev[1][0]
        elif kind=='release':
            for q in ev[1][1]:owner.pop(q,None)
    cur=None;pp={};lv=0
    for ev in b.wire_log:
        if ev[0]=='mark':cur=ev[1] if len(ev)>1 else cur
        elif ev[0]=='alloc':lv=ev[3];pp[str(cur)]=max(pp.get(str(cur),0),lv)
    json.dump(pp,open(HERE/f'parameters/{output}_phasepeaks.json','w'))
    print('wire_log kinds',dict(kinds));print('PEAK SNAPSHOT',dict(snap) if snap else None,flush=True)
    hook=dict(index=len(b._items),sigma=[list(w) for w in bm.SIG],tail=[list(w) for w in reversed(syn[35:])],
              free=list(b.free_qubits),active=bm.active,stopidx=list(bm.stopidx),parkidx=list(bm.parkidx))
    bm.reverse()
    assert not b._items._objs
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    np.save(HERE/f'parameters/{output}_tape.npy',rows)
    hook.update(maxq=b.next_qubit,peak=b.peak_qubits,peak_phase=b.peak_phase,outcomes=b.next_out,bits=b.next_bit,ops=len(rows),toffoli=b._stream_totals()[0],phases=marks)
    (HERE/f'parameters/{output}_hook.json').write_text(json.dumps(hook,indent=2))
    print({k:v for k,v in hook.items() if k not in ['sigma','tail','free','stopidx','parkidx','phases']},'seconds',time.time()-start,flush=True)
    json.dump(PEAKS,open(HERE/f'parameters/{output}_peaks.json','w'))
    TOFS.append(('end',b._stream_totals()[0]));json.dump(TOFS,open(HERE/f'parameters/{output}_tofs.json','w'))

if __name__=='__main__':
    import sys
    capture(sys.argv[1] if len(sys.argv)>1 else 'stream')
