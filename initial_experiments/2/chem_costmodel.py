from dfthcsos.verify import heuristic_costs
def est(R,B,C,lam=48.1,brot=12,J=9):
    L=76+R*B; K=-(-75//J)
    ang=2*K*L*1.012
    h=heuristic_costs(76,R,B,C,brot)
    T=ang+h['toffoli']; Q=h['persistent_qubits']+J*brot+max(h['temporary_qubits'],8)
    return lam*T*Q,T,Q,h
base=2903278577.4
for R,B,C in [(65,60,29),(65,40,29),(50,40,29),(40,40,29),(65,30,29),(40,30,29),(30,30,29),(50,24,29),(65,60,20),(40,40,40),(30,40,40)]:
    s,T,Q,h=est(R,B,C);print(R,B,C,'score %.3g'%s,'gain %.2f'%(base/s),int(T),Q,h['toffoli'],h['temporary_qubits'])
