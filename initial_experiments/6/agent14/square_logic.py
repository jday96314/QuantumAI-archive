"""Generate rounded-square truth tables; inputs 257..511 are don't-cares."""
from pathlib import Path
n=9;w=11
rows=['.i 9','.o 11','.ilb '+' '.join(f'a{i}' for i in range(n)),'.ob '+' '.join(f'o{i}' for i in range(w))]
for u in range(257):
 h=(u*u+32)//64
 rows.append(''.join(str(u>>j&1) for j in range(n))+' '+''.join(str(h>>j&1) for j in range(w)))
rows+=['.e']
Path('parameters/square.pla').write_text('\n'.join(rows)+'\n')
