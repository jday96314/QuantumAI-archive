"""Complete circuit scheduling search for the exact Karatsuba row rewrite."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import affine_invsqrt as affine
import binary_centered_invsqrt as poly
import fixed_top_affine_invsqrt as fixed
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact

HERE = Path(__file__).resolve().parent
SYNTHESIS_WORK = 250.85021091069865


def build(limit=90, keep=4, delay=2, product=17, deferred=3):
    old_shift, old_restore = poly.cshift, affine.restore_affine
    poly.cshift, affine.restore_affine = zero.cshift, fixed.restore_affine
    try:
        return affine.build(limit, keep, delay, product, deferred)
    finally:
        poly.cshift, affine.restore_affine = old_shift, old_restore


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wide', action='store_true')
    args = parser.parse_args()
    best = 182956.5271606425
    rows = []
    start = time.time()
    configurations = [(limit, keep, delay, product, deferred)
                      for limit in ([86,88,90,92,94,96,98,100] if args.wide else [88,90,92,94])
                      for keep in [2,3,4,5,6]
                      for delay in [2,3,4]
                      for product in [12,17,None]
                      for deferred in [2,3,4]]
    with (HERE/'results/row_search.jsonl').open('w') as log:
        for limit, keep, delay, product, deferred in configurations:
            parameters = dict(limit=limit, keep=keep, delay=delay, product=product, deferred=deferred)
            try:
                circuit = build(**parameters)
                circuit = compact(circuit, HERE/'candidates/row_temp.kmb')
            except (ValueError, AssertionError) as error:
                row = dict(parameters=parameters,status='INFEASIBLE',error=str(error))
                log.write(json.dumps(row)+'\n'); log.flush()
                continue
            tof = circuit.operation_statistics().num_toffoli_like
            cost = circuit.num_qubits*(tof-11+SYNTHESIS_WORK)
            row = dict(parameters=parameters, Q=circuit.num_qubits, static_Toffolis=tof,
                       projected_cost=cost, projected_overall=9.811059643318059*(182956.5271606425/cost)**(1/3))
            log.write(json.dumps(row)+'\n'); log.flush(); rows.append(row)
            if cost < best:
                best = cost
                path=HERE/'candidates/rows_best.kmb'
                circuit.to_file(str(path), format='kmb')
                row['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
                (HERE/'results/rows_best_build.json').write_text(json.dumps(row,indent=2)+'\n')
                print(json.dumps(dict(**row,seconds=time.time()-start)),flush=True)
    rows.sort(key=lambda row:row['projected_cost'])
    (HERE/'results/row_frontier.json').write_text(json.dumps(rows[:30],indent=2)+'\n')
    print(json.dumps(dict(status='SEARCH_COMPLETE',configurations=len(configurations),
                          feasible=len(rows),seconds=time.time()-start,best=rows[0] if rows else None)),flush=True)


if __name__ == '__main__':
    main()
