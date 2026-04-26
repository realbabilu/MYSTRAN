# Banded Optimization Proposal V1

## Purpose

This note proposes a staged optimization plan for the **banded solver subsystem only** in MYSTRAN.

The goal is to improve:
- memory efficiency
- banded solve performance
- conversion overhead from the default matrix path into banded storage
- solver robustness for medium-to-large bridge/frame/shell models

without disrupting the existing sparse solver ecosystem.

## Hard Scope Boundary

This work is **banded-only**.

Out of scope:
- changing the existing `CSR/CSC` sparse storage design
- changing the sparse assembly path used for `SuperLU`
- changing the `SuperLU_MT` path
- forcing sparse and banded solvers to share one internal matrix format
- replacing the sparse solver default behavior for large models

Allowed interaction with the sparse path:
- read existing matrix data that already exists in default storage
- make a dispatch decision about whether a problem is better treated as banded or sparse
- reuse ordering metadata such as `RCM` results when helpful

The sparse path remains:
- `CSR/CSC`
- `SuperLU`
- `SuperLU_MT`

unchanged in behavior and storage model.

## Current Situation

Today the banded path benefits from:
- the new `RCM` optimization
- bypassing legacy `BANDIT` in the preferred path
- existing LAPACK banded routines already present in the codebase
- partial SPD-aware logic already introduced in the LAPACK banded flow

However, the banded subsystem still has structural problems:
- the default matrix world is not banded-first
- conversion into band storage can be expensive
- conversion may involve extra transpose/reformat steps because the native matrix path is sparse-oriented
- fixed-width band storage can explode in memory if the bandwidth is still large after resequencing
- the code may rebuild banded forms more often than necessary
- SPD opportunities are not yet treated as a first-class optimized path

## Main Design Principle

We should treat banded as a **specialized acceleration path**, not as the universal matrix core.

That means:
1. keep sparse general and untouched
2. make banded smarter about when it is selected
3. make banded cheaper when selected
4. avoid premature introduction of new storage formats until the current banded path is fully understood

## Optimization Goals

### Primary goals
- Reduce peak RAM for banded solves.
- Reduce banded conversion cost from the default matrix path.
- Reduce repeated format conversions.
- Reuse banded metadata and workspaces when possible.
- Exploit SPD cases safely with half-band storage.

### Secondary goals
- Improve portability across:
  - Windows PC
  - Linux PC
  - ARM64
  - Apple Silicon
- Avoid dependence on vendor-specific math stacks such as MKL.
- Keep future compatibility with optional external solvers such as PARDISO, without designing around them now.

## Non-Goals

This proposal does **not** attempt to:
- make banded beat sparse on every model
- replace sparse for huge irregular models
- force skyline/profile storage immediately
- rewrite the global assembly engine in one step

## Proposed Phases

## Phase 0: Audit and Instrumentation

### Objective
Understand the current banded pipeline precisely before making structural changes.

### Questions to answer
- Where does MYSTRAN decide that a matrix will use a banded path?
- Where does conversion from default storage to band storage happen?
- How many times is the same banded matrix built or rebuilt in one run?
- Where are transposes or symmetric-to-general reshapes happening?
- What are the largest allocations in the banded path?
- Which parts are conversion cost, factor cost, and solve cost?

### Deliverables
- a call-flow note for the current banded path
- timing instrumentation for:
  - ordering
  - convert-to-band
  - factor
  - solve
- memory estimates printed to `.F06` or debug stream:
  - full-band estimate
  - half-band estimate
  - possible skyline estimate later

### Why this phase matters
Without this, we risk optimizing the wrong layer.

## Phase 1: Banded Feasibility and Dispatch

### Objective
Make MYSTRAN smarter about when banded is worthwhile.

### Proposed work
- Keep `RCM` as the preferred banded ordering path.
- Add bandwidth reporting before and after `RCM`.
- Add estimated storage reporting:
  - general band
  - SPD half-band candidate
- Add a banded feasibility check:
  - if post-RCM bandwidth is still too large, prefer sparse
  - if the model shape is favorable, allow banded

### Expected benefit
- Prevent catastrophic banded RAM growth.
- Avoid forcing banded on models that are structurally sparse but not truly narrow-band.

## Phase 2: SPD Half-Band Optimization

### Objective
Exploit safe SPD cases for immediate memory savings and faster factorization.

### Proposed work
- Introduce a clear `SPD candidate` classification for the banded path.
- Use half-band storage for validated SPD banded problems.
- Route that path to:
  - `DPBTRF`
  - `DPBTRS`
- If the SPD assumption fails during factorization, fall back safely to:
  - general banded, or
  - sparse path

### Important caution
Do **not** assume all linear static/modal problems are automatically SPD.

Possible SPD breakers include:
- penalties
- releases
- constraints
- drilling stabilizations
- formulation-specific terms

So the logic should be:
- `candidate first`
- `validated during factorization`
- `fallback if rejected`

### Expected benefit
- About 50 percent band-storage savings in good SPD cases
- Faster factor/solve than generic band LU

## Phase 3: Conversion Path Optimization

### Objective
Reduce the overhead of building banded matrices from the default storage world.

### Proposed work
- Identify all current conversion and transpose steps.
- Remove redundant reformatting where possible.
- Introduce a compact banded metadata object containing:
  - bandwidth
  - storage class (`DGB` or `DPB`)
  - symmetry class
  - SPD candidate flag
  - row/column layout expectations
- Build banded matrices in a single planned conversion whenever possible.

### Key rule
Do not rewrite the sparse matrix representation.
Instead, optimize the bridge from sparse-oriented storage into banded storage.

### Expected benefit
- Lower conversion penalty
- Cleaner code paths
- Better basis for later direct-to-band assembly

## Phase 4: Reuse of Banded Work and Factors

### Objective
Stop rebuilding expensive banded artifacts unnecessarily.

### Proposed work
- Reuse:
  - banded descriptors
  - work arrays
  - factor workspace sizing
  - ordering metadata
- For repeated solves, avoid full rebuild when matrix structure is unchanged.
- For modal workflows, explicitly support reuse where mathematically valid.

### Expected benefit
- Lower total runtime in multi-step or iterative workflows
- Better user experience for repeated analysis cases

## Phase 5: Direct-to-Band Assembly Feasibility

### Objective
Evaluate whether element assembly can target banded storage directly after final DOF ordering is known.

### Why this matters
This is the first major step that could remove a large part of conversion overhead.

### Proposed work
- Study whether element contributions can be assembled directly into:
  - general band arrays
  - half-band SPD arrays
- Keep this isolated from the sparse path.
- If direct-to-band assembly is added, it should be enabled only for the banded branch after dispatch selection.

### Risks
- More invasive than Phases 1 to 4
- Harder to maintain
- Must not leak assumptions into sparse assembly

### Expected benefit
- Remove convert/re-convert cost almost entirely for true banded jobs

## Phase 6: Skyline/Profile Evaluation

### Objective
Decide with evidence whether skyline/profile storage is worth implementing.

### Recommendation
Do **not** assume skyline is automatically the next step.

Evaluate first:
- full-band bytes
- half-band bytes
- skyline/profile bytes
- conversion cost
- factor cost
- implementation complexity
- maintainability cost

### Decision rule
Only move to skyline if real bridge/shell benchmarks show a substantial advantage over:
- RCM + half-band
- optimized general band
- direct-to-band assembly

### Why delay skyline
- It is a large architectural branch.
- LAPACK support is less straightforward than banded.
- It may add long-term maintenance cost for limited real gain.

## Phase 7: High-Performance Banded Kernel Work

### Objective
Optimize the final chosen banded storage for modern CPUs.

### Possible work
- cache-aware blocking
- reduced-copy kernels
- active-column or frontal-style local optimizations
- OpenMP where safe
- ARM64- and Apple-Silicon-friendly tuning

### Constraint
This phase should happen only after the storage strategy is stable.

## Proposed Metrics

For each benchmark case, collect:
- matrix size
- original bandwidth
- post-RCM bandwidth
- general-band bytes
- half-band bytes
- skyline-estimated bytes if available
- convert time
- factor time
- solve time
- repeated-solve time
- total runtime

## Suggested Benchmarks

Use a mix of:
- bridge-like beam/frame models
- shell models
- medium symmetric models likely to favor banded
- models that clearly should remain sparse

At minimum:
- one Eigen benchmark already used in current work
- one or more static plate/shell cases
- one bridge-style frame case with significant bandwidth pressure

## Recommended Immediate Next Steps

1. Audit the current banded path end-to-end.
2. Add instrumentation for bandwidth, storage estimates, and convert/factor/solve timing.
3. Implement a clean SPD-candidate path for banded only.
4. Prototype half-band storage using existing `DPBTRF/DPBTRS`.
5. Measure again before considering skyline.

## Summary

The best near-term strategy is not a full storage rewrite.

The best near-term strategy is:
- keep sparse untouched
- optimize banded selection
- optimize banded conversion
- exploit SPD half-band aggressively but safely
- reuse work/factors
- postpone skyline until measurements justify it

That path is lower risk, portable, and much more likely to deliver practical wins quickly.
