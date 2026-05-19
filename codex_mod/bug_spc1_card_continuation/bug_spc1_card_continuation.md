# SPC1 Continuation Truncation

This package captures the `SPC1` card-continuation bug and its local source
hardening in the current MYSTRAN tree.

## Problem

The generated benchmark decks originally wrote too many grid IDs on a single `SPC1` line, for example:

```text
SPC1,1,123456,1,2,3,4,5,6,7,8,9
```

That format is not valid as an unlimited one-line free-field list in MYSTRAN.

The bulk-data reader only accepts:

- parent card: grid IDs in fields `4-9`
- continuation cards: grid IDs in fields `2-9`

Source proof:

- [BD_SPC1.f90](C:/mystran3/MYSTRANSolver-18.0.0/Source/LK1/L1A-BD/BD_SPC1.f90)
- [FFIELD.f90](C:/mystran3/MYSTRANSolver-18.0.0/Source/LK1/L1A/FFIELD.f90)

Current snippet markers used for this issue:

```fortran
! --- spc1_fix begin --- !
! --- spc1_fix end --- !
```

See:

- [bd_spc1_continuation_snippet.f90](C:/mystran3/codex_mod/bug_spc1_card_continuation/bd_spc1_continuation_snippet.f90)
- [BD_SPC1.f90](C:/mystran3/codex_mod/bug_spc1_card_continuation/BD_SPC1.f90)

## Wrong vs Correct Deck Form

Wrong single-line form:

```text
SPC1,1,123456,1,2,3,4,5,6,7,8,9
```

Correct continued form:

```text
SPC1,1,123456,1,2,3,4,5,6
,7,8,9
```

or, more compactly:

```text
SPC1,1,123456,1,THRU,9
```

## Local hardening behavior

The local source now rejects the malformed one-line case explicitly instead of
silently truncating it:

- `BD_SPC1` raises `*ERROR 1129` when trailing inline data reaches field 10
- `FFIELD` raises `*ERROR 1003` when a free-field physical line exceeds the
  supported field count

This means the old hidden partial-clamp behavior is now converted into a
visible input error.

## Proof It Truncated

The strongest proof is behavioral:

- before fixing continuation, base grids `7-9` showed nonzero displacement in the MYSTRAN static unit-load run
- after fixing continuation, all fixed-base grids returned exactly zero displacement
- after the fix, MYSTRAN matched OpenSees node-for-node for the full `1-024` static displacement profile

## Numerical Evidence

Before the `SPC1` fix:

- node `29` unit-load ratio was `1.321610`
- node `19` unit-load ratio was `1.190371`
- floor `z=13` average displacement ratio was `1.590347`
- base grids `7-9` were moving even though they should have been fixed

After the `SPC1` fix:

- node `29`: OpenSees `1.714145015e-04`, MYSTRAN `1.714145000e-04`
- node `19`: OpenSees `1.645314940e-04`, MYSTRAN `1.645315000e-04`
- story `z=13` average ratio `1.000000`
- story `z=26` average ratio `1.000000`
- modal periods: `0.227062, 0.215633, 0.073345, 0.072005 s`

## Conclusion

The `SPC1` continuation issue is real and now closed locally:

- malformed long one-line `SPC1` free-field input is rejected
- valid continuation-form `SPC1` still works
- valid `THRU` form still works
- this prevents silent partial base-fixity loss in downstream benchmarks

## Included Files

- [bd_spc1_continuation_snippet.f90](C:/mystran3/codex_mod/bug_spc1_card_continuation/bd_spc1_continuation_snippet.f90)
- [bad_vs_good_spc1_cards.txt](C:/mystran3/codex_mod/bug_spc1_card_continuation/bad_vs_good_spc1_cards.txt)
- [BD_SPC1.f90](C:/mystran3/codex_mod/bug_spc1_card_continuation/BD_SPC1.f90)
