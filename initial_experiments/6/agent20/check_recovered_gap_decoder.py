"""Native restoration/phase checks for the locally recovered gap circuit."""
from pathlib import Path
import hashlib
import json
import numpy as np
import recovered_gap_decoder as research

HERE = Path(__file__).resolve().parent
km = research.km


def main():
    n = 14
    path = HERE/f'candidates/recovered_gap_n{n}_both.kmb'
    circuit = km.Circuit.from_file(str(path))
    f = km.GF2Field(12)
    rng = np.random.default_rng(2026102531)
    samples, labels = [], []
    for weight in range(8):
        for _ in range(64):
            roots = rng.choice(np.arange(1,4096), weight, replace=False)
            magnitudes = rng.integers(1,4096,weight)
            syndrome = [0]*n
            for x,e in zip(roots,magnitudes):
                power = int(x)
                for i in range(n):
                    syndrome[i] ^= f.mul(int(e), power)
                    power = f.mul(power, int(x))
            samples.append(syndrome)
            labels.append(f'physical_weight{weight}')
    for _ in range(512):
        samples.append(rng.integers(0,4096,n).tolist())
        labels.append('uniform_syndrome')
    # Deliberately unsupported zero patterns are reported separately.
    for i in range(256):
        sample = rng.integers(0,4096,n).tolist()
        sample[0] = 0
        if i%2:
            sample[1] = 0
        samples.append(sample)
        labels.append('forced_zero_chart')
    data = np.array(samples, dtype=np.uint64)
    registers = [circuit.register_data[f'S{i}'] for i in range(n)]
    public = {q for w in registers for q in w}
    scratch = km.array([km.q(i) for i in range(circuit.num_qubits) if km.q(i) not in public])
    bit = lambda a: ((a[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    failures = []
    for seed in [2026102531,2026102532]:
        for start in range(0,len(data),256):
            values = data[start:start+256]
            sim = km.Simulator(batch_size=len(values), seed=seed+start)
            sim.use_same_registers_as(circuit)
            for i,w in enumerate(registers):
                sim.write_across_shots(w,bit(values[:,i]))
            sim.do(circuit)
            value_bad = np.zeros(len(values),dtype=bool)
            for i,w in enumerate(registers):
                value_bad |= (sim.read_across_shots(w)!=bit(values[:,i])).any(0)
            scratch_bad = sim.read_across_shots(scratch).any(0)
            phase_bad = sim.read_phase_flipped_across_shots()
            for k in np.flatnonzero(value_bad|scratch_bad|phase_bad):
                failures.append(dict(seed=seed,index=start+int(k),label=labels[start+k],value=bool(value_bad[k]),scratch=bool(scratch_bad[k]),phase=bool(phase_bad[k])))
    report = dict(status='DIAGNOSTIC_COMPLETE',Q=circuit.num_qubits,T=circuit.operation_statistics().num_toffoli_like,
                  shots=2*len(data),failures=failures,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  scope='paired restoration and measurement-phase only; no Forney truth check or full benchmark acceptance')
    (HERE/'results/recovered_gap_n14_native_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
