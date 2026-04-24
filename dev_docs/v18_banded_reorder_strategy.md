# v18 Banded Reorder Strategy (Add-on vs Classic BANDIT)

Date: 2026-04-22  
Status: Proposed for implementation

! !--- RCM BANDED ADD-ON --- begin!
## Context

On `midas_static_24_cquad4_moremesh` we compared three `SOLLIB=BANDED` paths:

1. Baseline banded (`INPUT` order)
2. Classic BANDIT (`PARAM,GRIDSEQ,BANDIT`)
3. Add-on RCM reorder + sorted GRID output (external pre-process)

## Measured Runtime (6 runs average)

| Path | Avg time |
|---|---:|
| Baseline banded | 2.6109 s |
| Classic BANDIT + banded | 2.6852 s |
| Add-on RCM + sorted + banded | 0.4696 s |

Observation: add-on RCM is about 5.5x faster than baseline and clearly faster than classic BANDIT for this case.

## Decision

Use **add-on reorder strategy** as primary optimization for banded solver.  
Do **not** replace/remove BANDIT in v18 now.

Rationale:
- Lowest risk to core solver behavior.
- Keeps existing BANDIT path for backward compatibility and fallback.
- Add-on approach already proven by benchmark and can be introduced behind a parameter gate.

## Proposed Integration (In-Core)

### New parameter

- `PARAM,BANDEDOPT,YES/NO` (default `NO`)
- Active only when `SOLLIB=BANDED`

### Execution flow

1. Parse input as usual.
2. If `SOLLIB=BANDED` and `BANDEDOPT=YES`:
   - Build graph from structural connectivity.
   - Run RCM node ordering.
   - Rebuild internal GRID sequence before DOF table finalization.
3. Continue existing banded assembly/factorization path.
4. If optimization step fails, print warning and fall back to existing flow.

## Scope and Non-Goals

In scope:
- Banded runtime reduction via ordering only.
- No element formulation changes.
- No sparse solver behavior changes.

Out of scope (this phase):
- Replacing BANDIT implementation.
- Automatic SPD detection/switching in sparse path.
- Eigen04 mechanism fixes via solver trick (those remain boundary-condition/modeling tasks).

## Rollout Plan

1. **Phase A**: parameter plumbing (`PARAMS`, `BD_PARAM`) + diagnostics print.
2. **Phase B**: internal RCM ordering hook at Link0 pre-DOF-final stage.
3. **Phase C**: regression run:
   - static/eigen existing validation decks
   - static-24 performance comparison (baseline vs `BANDEDOPT=YES`)
4. **Phase D**: document results in validation report and keep BANDIT fallback.

## Acceptance Criteria

- Numerical outputs match baseline (within normal floating tolerance).
- No increase in fatal errors across existing validation suite.
- Runtime improvement on static-24 style mesh case is significant (target >= 2x).

## Risk Notes

- Reordering must preserve all references (loads, SPC, MPC, recovery indexing).
- Any failure in ordering path must be non-fatal and fall back automatically.
- Keep output transparency: report old/new bandwidth proxy and elapsed time.
! !--- RCM BANDED ADD-ON --- end!
