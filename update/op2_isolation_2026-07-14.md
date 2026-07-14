# OP2 Isolation Notes - 2026-07-14

## 1. `vic/7/S30 mitc4 thermal strain free diamond cylinder.bdf`

### F06 path

- Current F06-based validation for this deck passes on the active checks.
- The previously reported residual `Error = 2.5` does not reproduce on the current F06 path.

### OP2 path

- Direct OP2 probing shows that this deck is **not yet safe as an OP2 parity case** for shell result validation.
- The reference MSC OP2 and the MYSTRAN OP2 both load through `pyNastran`, but several values differ strongly from the F06/cases semantics:

Reference MSC OP2:

- `SC/1/SHELLFORCES/EID/174/CORNER/2/QX = 144169.2`
- `SC/1/SHELLFORCES/EID/174/CORNER/2/QY = 144169.2`
- `SC/1/SHELLSTRESSES/EID/190/CORNER/1/Z2/VONMISES = 24564700.0`
- `SC/1/SHELLSTRESSES/EID/243/CORNER/3/Z2/XY = -2646352.0`

Isolated MYSTRAN OP2:

- `SC/1/SHELLFORCES/EID/174/CORNER/2/QX = 20392.63`
- `SC/1/SHELLFORCES/EID/174/CORNER/2/QY = -95758.51`
- `SC/1/SHELLSTRESSES/EID/190/CORNER/1/Z2/VONMISES = 96104140.0`
- `SC/1/SHELLSTRESSES/EID/243/CORNER/3/Z2/XY = 1301795.0`

### Interpretation

- This is not just a MYSTRAN writer issue.
- `cases.txt` values for this deck were authored around F06 semantics.
- The MSC reference OP2 shell-force/shell-stress values do not match those same `cases.txt` values either.
- Therefore this deck should **not** be used as an OP2 parity oracle until the shell-force/shell-stress OP2 semantics are explicitly aligned and documented.

### Additional gaps

- `SC/1/SHELLSTRAINS/GID/350/Z1/YY` is `None` from OP2 while available through F06 node averaging.
- `SC/1/SHELLSTRESSES/EID/322/CORNER/3/ZMID/VONMISES` is also absent in the current OP2 path.
- These are expected: OP2 gives element/corner payloads, while some F06 validations rely on node-averaged or derived views.


## 2. MSC RSA reference: `femap_RSA/msc/beam_rs_msc2.bdf.op2`

### Input

- Deck file: `D:\18a\femap_RSA\msc\beam_RS_msc2.bdf.txt`
- OP2 file: `D:\18a\femap_RSA\msc\beam_rs_msc2.bdf.op2`

### Deck characteristics

- `SOL 103`
- `DISPLACEMENT(PLOT)=ALL`
- `VELOCITY(PLOT)=ALL`
- `ACCELERATION(PLOT)=ALL`
- `SPCFORCE(PLOT)=ALL`
- `STRESS(PLOT,CORNER)=ALL`
- `MEFFMASS(ALL)=YES`
- `PARAM,SCRSPEC,0`
- `PARAM,OPTION,SRSS`

### Parser blockers found

Before patching `op2_query.py`, pyNastran aborted on:

1. `OUPV1` velocity with `thermal=4`
2. `OUG` acceleration with `thermal=4`
3. `OQP1/OQG` SPCFORCE with `thermal=4`

### Defensive parser patches added

In `D:\18a\MYSTRAN_Validation-main\op2_query.py`:

- skip unsupported OUG velocity subtable when `thermal=4`
- skip unsupported OUG acceleration subtable when `thermal=4`
- skip unsupported OQG/SPCFORCE subtable when `thermal=4`

These patches are intentionally non-destructive:

- they do **not** invent RSA semantics
- they only prevent one unsupported RSA subtable from aborting the entire OP2 read
- other tables can continue loading

### Current readback status after patch

- `SC/1/REALEIGENVALUES/MODE/1/EIGENVALUE = 1017.487`
- RSA displacement/velocity/acceleration/SPCFORCE are still `None` in `OP2Query`

### Interpretation

- We now have a stable partial-open path for MSC RSA OP2.
- The next step is not more skip patches; it is to implement actual `thermal=4` RSA result decoding/mapping for:
  - displacement
  - velocity
  - acceleration
  - SPC forces

Until that is done, `beam_rs_msc2.bdf.op2` is useful for:

- table presence audit
- eigenvalue audit
- confirming which RSA result families still need dedicated mapping

but **not yet** for full scalar-response validation through `op2_query.py`.
