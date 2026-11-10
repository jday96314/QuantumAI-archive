# Round 5: Five targeted Codex experiments

I would run **three agents initially**, with two backup experiments if any of the first three reach an unpromising result.

The ranking emphasizes the probability of achieving a substantial, **verified** improvement within approximately one round-4-sized agent effort, rather than maximum theoretical upside.

| Rank | Working directory | Experiment | Priority |
|---|---|---|---|
| **1** | `5/dqi-repair` | Repair exceptional decoding failures in the compact DQI decoder | Run immediately |
| **2** | `5/chem-qfit` | Develop a better quantization-aware optimizer and compress the chemistry factorization | Run immediately |
| **3** | `5/dqi-phase` | Replace expensive DQI phase extraction with a width-efficient, division-free method | Run immediately |
| 4 | `5/dqi-hybrid` | Combine compact Euclidean storage with cheaper Berlekamp–Massey arithmetic | Backup |
| 5 | `5/chem-blocks` | Explore resource-optimized multi-block chemistry factorizations | Backup |

**These are deliberately separate experiments**, although their successful components could eventually be combined.

I rank DQI exception repair first because it has a relatively clear mathematical target and an approximately 1.42× theoretical improvement available from correctness alone. The chemistry experiments have greater potential to escape the existing architecture, but considerably more numerical-optimization risk.

For reference, the starting point is the verified round-4 result:

| Task | Current raw cost |
|---|---:|
| Inverse square root | 182,956.5272 |
| DQI | 1,598,204,229,385 |
| Chemistry | 708,547,388.96 |
| **Combined score** | **9.367722** |

A task-level improvement of 1.5× increases the combined score by approximately 14.5%, assuming the other tasks remain unchanged.

The prompts below distinguish complete verifier-confirmed improvements from component-level results and resource projections. They also require numerical controls so that a failed optimizer isn't mistaken for proof that an architecture is infeasible.

## Ready-to-use Codex prompts

I've written five standalone prompts. Each includes the technical hypothesis, baseline measurements, implementation strategy, validation requirements, success criteria, and conditions for abandoning an unproductive approach.

**:chatgpt-content-reference{index="2"}[Download all five Codex prompts in one Markdown file](sandbox:/mnt/data/round5_codex_prompts.md)**

### Initial three agents

1. **`dqi-repair`** — round5_dqi_repair_prompt.md[Codex prompt](sandbox:/mnt/data/round5_dqi_repair_prompt.md)  
   Diagnose the 380 phase errors, identify their underlying cause, and implement reversible exception handling. This is the most immediately actionable opportunity.

2. **`chem-qfit`** — round5_chem_qfit_prompt.md[Codex prompt](sandbox:/mnt/data/round5_chem_qfit_prompt.md)  
   Reproduce the accepted quantized incumbent, develop better discrete-angle optimization, and use it to make smaller chemistry factorizations viable.

3. **`dqi-phase`** — round5_dqi_phase_prompt.md[Codex prompt](sandbox:/mnt/data/round5_dqi_phase_prompt.md)  
   Investigate division-free phase extraction, especially the in-place transformation \(F=\omega\Psi\bmod\sigma\), while accounting for the complete reversible circuit and its peak qubit count.

### Two backup agents

4. **`dqi-hybrid`** — round5_dqi_hybrid_prompt.md[Codex prompt](sandbox:/mnt/data/round5_dqi_hybrid_prompt.md)  
   Combine the low gate counts of the Berlekamp–Massey decoder with the substantially smaller workspace achieved in round 4.

5. **`chem-blocks`** — round5_chem_blocks_prompt.md[Codex prompt](sandbox:/mnt/data/round5_chem_blocks_prompt.md)  
   Search the exact chemistry resource frontier, then investigate physically valid multi-block factorizations, including candidates that exploit seven-bit QROAM width discontinuities.

## Additional research directions

**:chatgpt-content-reference{index="8"}[Download: Round 5 — Other Promising Directions](sandbox:/mnt/data/round5_other_directions.md)**

This document covers ten additional opportunities, including content-aware ANGLE compilation, completing the Wronskian decoder, inverse-square-root arithmetic redesign, Givens-tree optimization, and exploiting the classical DQI masks.

For each, it explains the proposed experiment, what would constitute meaningful evidence, and why I ranked it below the initial five.

## How I'd manage the round

All five prompts instruct Codex to preserve the existing verified submission, work independently, and maintain a `RESULTS.md` containing measured results, failed experiments, and reproduction instructions.

There are two especially important safeguards.

**First, DQI failure reduction must be evaluated using the complete cost formula.** An exception-handling circuit that eliminates almost all failures can still lose if its reversible bookkeeping increases width or gate count too much. The current 961-qubit decoder and its 380 failures provide a concrete starting point. REPORT

**Second, chemistry experiments must distinguish continuous fitting, integer-angle fitting, and complete native verification.** Round 4 demonstrated that the quantization pipeline itself could fail to reproduce a previously successful solution. REPORT

The prompts also prohibit agents from examining other round-5 workspaces, while allowing them to reuse prior-round implementations and public research.

### One scheduling recommendation

If `dqi-repair` succeeds, I would **not immediately merge its solution into `dqi-phase` or `dqi-hybrid` while those agents are running**.

Instead, preserve the independent experiments. Once the agents finish, integrate the successful changes into a separate candidate and verify their combined cost.

That preserves your ability to determine whether each idea actually works and avoids disrupting an agent's partially completed implementation.

**The three initial experiments provide a useful balance:** one relatively focused correctness improvement, one numerical algorithm breakthrough attempt, and one substantially different DQI phase-extraction architecture.