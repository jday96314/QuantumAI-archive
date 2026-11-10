"""Unmodified host PREPARE plus explicit generalized signed ANGLE bounds.

Reuses completed chem-blocks/frontier.py price, whose ANGLE estimates are
constructive worst-case bounds (not a physical or native submission score).
"""
import csv
from core import *
from frontier import price

def main():
    rows=[];ev=native_spectrum();positive=ev[ev>0];negative=float(np.square(ev[ev<0]).sum())
    for R,B in [(4,180),(4,224),(2,351),(2,375),(2,384),(2,447),(2,479),(3,240),(3,272),(4,256),(6,128)]:
      for C in sorted(set([88,96,112,128,144,152,160,168,175,180,184,188,189,190,191,192,196,200,208,224,256,300,350,B])):
       if C>B:continue
       for bits in [8,9]:
        options=[price(R,B,C,bits,J) for J in range(4,17)]
        xs=[x for x in options if x]
        if not xs:continue
        x=min(xs,key=lambda x:x['conditional_cost_at_incumbent_lambda'])
        x['lambda_max_for_1p3']=708547388.9553173/(1.3*x['T_total_bound']*x['Q_total_bound'])
        x['lambda_max_for_1p5']=708547388.9553173/(1.5*x['T_total_bound']*x['Q_total_bound'])
        x['unrestricted_PSD_rank_error_lower_bound']=float(np.sqrt(negative+np.square(positive[R*C:]).sum()))
        x['necessary_rank_test_pass']=x['unrestricted_PSD_rank_error_lower_bound']<=.02
        if not x['necessary_rank_test_pass']:x['physical_feasibility']='EXCLUDED_BY_UNRESTRICTED_PSD_RANK'
        rows.append(x)
    with (HERE/'results/resources.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    report(HERE/'results/resources.json',dict(status='CONDITIONAL_RESOURCE_BOUNDS_ONLY',rows=rows))
if __name__=='__main__':main()
