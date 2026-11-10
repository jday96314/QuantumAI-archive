"""Aggregate actual allocator lifetimes without retaining a massive event log."""
from collections import Counter


def instrument(A):
    opened=A.open;closed=A.close;live={};summary={};peaks=[];b=A.b
    def start(name,n):
        before=b.peak_qubits;word=opened(name,n);live[tuple(word)]=(name,n)
        row=summary.setdefault(name,dict(opens=0,max_width=0,peak_live=0));row['opens']+=1;row['max_width']=max(n,row['max_width']);row['peak_live']=max(b.active_qubits,row['peak_live'])
        if b.peak_qubits>before:
            peaks.append(dict(peak=b.peak_qubits,cause=name,windows=[dict(name=nm,width=width) for nm,width in live.values()]))
        return word
    def end(name,word):
        closed(name,word);live.pop(tuple(word),None)
    A.open=start;A.close=end
    return lambda:dict(peak=b.peak_qubits,wire_high_water=b.next_qubit,summary=summary,new_peaks=peaks,remaining_windows=[dict(name=nm,width=width) for nm,width in live.values()])
