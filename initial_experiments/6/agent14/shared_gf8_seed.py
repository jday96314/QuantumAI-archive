"""Quadratic seed for two GF8 products with a common first operand."""
from pathlib import Path
import kickmix as km
from crt_multiplier import small_terms
terms=small_terms(km.GF2Field(3).modulus);out=[(a,b<<shift,c<<dest) for shift,dest in [(3,0),(6,3)] for a,b,c in terms];assert len(out)==12;Path('parameters/shared_gf8_seed.txt').write_text('12\n'+''.join(f'{a} {b} {c}\n' for a,b,c in out))
