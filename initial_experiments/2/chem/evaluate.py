import sys,json,numpy as np
from dfthcsos.verify import load_inputs,evaluate_words
z=np.load(sys.argv[1]);inp=load_inputs(sys.argv[1])
print(json.dumps(evaluate_words(inp,z['angle_words'],z['tree_edges'],12)))
