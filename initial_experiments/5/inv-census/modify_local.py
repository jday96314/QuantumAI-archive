"""Apply checked local-only experimental transformations to copied baseline."""
from pathlib import Path
H=Path(__file__).resolve().parent
p=H/'binary_centered_invsqrt.py';s=p.read_text()
def rep(a,b):
 global s
 assert s.count(a)==1,(a,s.count(a));s=s.replace(a,b)
rep('    def shift(inv):\n        if binary:', '    def shift(inv):\n        if globals().get("NEW_PHASE_MODE"):return\n        if binary:')
rep('values=(tab[:,0]>>cutoff)&((1<<delay_high_bits)-1);anc=y[sizes[0]:]',
    'values=(tab[:,0]>>cutoff)&((1<<delay_high_bits)-1);anc=y[sizes[0]:]\n        local_anc=b.alloc_qubits(max(0,len(address)-1-len(anc)))\n        anc=km.array([*anc,*local_anc])')
rep('visit(address[-1],len(address)-2,mid,len(values));b.x(x[20-frac])',
    'visit(address[-1],len(address)-2,mid,len(values));b.x(x[20-frac]);b.free(local_anc)')
rep('    else:nr._phase_by_bits(b,phase_y)',
    '''    else:
        mode=globals().get("NEW_PHASE_MODE")
        if mode=="borrowed":
            from exponent_kernels import borrowed_phase
            borrowed_phase(b,y[:27],x,flags)
        elif mode=="walsh":
            from exponent_kernels import walsh_phase
            from fractions import Fraction
            import json
            co={int(k):Fraction(v) for k,v in json.loads((Path(__file__).resolve().parent/'results/walsh_search.json').read_text())['walsh_coefficients'].items()}
            walsh_phase(b,y[:27],flags[1:],co)
        else:nr._phase_by_bits(b,phase_y)''')
p.write_text(s)
p=H/'affine_invsqrt.py';s=p.read_text()
rep('b.free(affine_extra);allocated_y=allocated_y[:-{deferred}];y=y[:-{deferred}]',
    'b.free(affine_extra);allocated_y=allocated_y[:-{deferred}] if {deferred} else allocated_y;y=y[:-{deferred}] if {deferred} else y')
rep('reuse_normalized_top=True,output_extra_bits=10,defer_workspace=True,',
    'reuse_normalized_top=True,output_extra_bits=0 if getattr(p,"NEW_PHASE_MODE",None) else 10,defer_workspace=True,')
p.write_text(s)
