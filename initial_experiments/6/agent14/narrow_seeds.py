from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
for w in [(27,17,8),(27,18,7),(26,18,8),(26,17,8)]:
 tab=np.rint(p.float_table(7,2)*np.exp2(w)).astype(np.int64)
 tag='_'.join(map(str,w));Path(f'parameters/narrow_seed_{tag}.txt').write_text(str(len(tab))+'\n'+''.join(' '.join(map(str,row))+'\n' for row in tab))
