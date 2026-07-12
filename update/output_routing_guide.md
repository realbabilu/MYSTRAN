# Output Routing Guide for `v18.00.a`

This note explains the current output-routing behavior in this branch for:

- `F06`
- `OP2`
- `NEU`

It focuses on one practical question:

- how to get `OP2` without detailed displacement/force/stress/strain tables in `F06`
- how to get `NEU` without those detailed tables in `F06`
- how to get all three
- how to get only one of them

## Short version

Use:

- `PARAM,OUTMODE,SMART`

when you want `F06` to behave more like a summary/diagnostic file and not a bulk-results dump.

Then use explicit Case Control descriptors such as:

- `DISP(PLOT)=ALL`
- `DISP(NEU)=ALL`
- `DISP(PRINT,PLOT,NEU)=ALL`

Do not rely on bare requests like:

- `DISP = ALL`
- `STRE = ALL`
- `STRN = ALL`
- `ELFO = ALL`
- `SPCF = ALL`
- `OLOAD = ALL`

because in the current branch these bare forms default to `PRINT+PLOT`, so they still send detailed result tables to `F06`.

## Main controls

Bulk Data `PARAM` flags:

- `PARAM,PRTF06,Y` forces broad `F06` output in `LEGACY`
- `PARAM,PRTOP2,Y` forces broad `OP2` output
- `PARAM,PRTNEU,Y` forces broad `NEU` output
- `PARAM,PRTALL,Y` turns on all three broad file targets
- `PARAM,OUTMODE,LEGACY` keeps old broad writer behavior
- `PARAM,OUTMODE,SMART` keeps diagnostics in `F06`, but avoids using `PRTF06` as a blanket “print all result tables” switch

Case Control descriptors per result family:

- `PRINT`
- `PLOT`
- `NEU`
- `CSV`
- `PUNCH` where supported

## What `OUTMODE=SMART` currently does

In this branch, `SMART` mainly affects the blanket `PRTF06` behavior in `LINK9`.

Practical effect:

- `PRTF06=Y` no longer auto-forces all result families into detailed `F06` tables
- `PRTOP2=Y` still enables broad `OP2` output
- `PRTNEU=Y` still enables broad `NEU` output
- `F06` still receives normal run text:
  - job progress
  - warnings
  - fatal/error text
  - summaries

Important boundary:

- some detailed `F06` result printing still depends on the Case Control family descriptors
- bare `... = ALL` requests are still expanded by several `CC_*` routines to `PRINT+PLOT`

So `SMART` helps, but it is not enough by itself if the deck still says `DISP = ALL`.

## Why bare `= ALL` is dangerous for `OP2-only`

Several Case Control processors in this branch explicitly convert a bare request into `PRINT+PLOT`.

Examples:

- [CC_DISP.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_DISP.f90)
- [CC_STRE.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_STRE.f90)
- [CC_STRN.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_STRN.f90)
- [CC_ELFO.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_ELFO.f90)
- [CC_SPCF.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_SPCF.f90)
- [CC_OLOA.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_OLOA.f90)

That means:

- `DISP = ALL` does not mean “let `OP2` decide”
- it means “print to `F06` and write to `OP2`”

If you want to avoid detailed `F06` result tables, use explicit descriptors.

## Recommended recipes

### 1. `OP2` only, with `F06` kept as summary/diagnostic text

Use:

```text
PARAM,OUTMODE,SMART
PARAM,PRTOP2,Y
PARAM,PRTF06,N
PARAM,PRTNEU,N
```

Case Control examples:

```text
DISP(PLOT)=ALL
SPCF(PLOT)=ALL
OLOAD(PLOT)=ALL
ELFO(PLOT)=ALL
STRE(PLOT)=ALL
STRN(PLOT)=ALL
```

Expected behavior:

- `OP2` gets the result payload
- `F06` remains as the normal report/summary/warning file
- detailed displacement/force/stress/strain tables should stay out of `F06` as long as you did not request `PRINT`

Do not use:

```text
DISP=ALL
STRE=ALL
STRN=ALL
```

because those bare forms default back to `PRINT+PLOT`.

### 2. `NEU` only, with `F06` kept as summary/diagnostic text

Use:

```text
PARAM,OUTMODE,SMART
PARAM,PRTOP2,N
PARAM,PRTF06,N
PARAM,PRTNEU,Y
```

Case Control examples:

```text
DISP(NEU)=ALL
SPCF(NEU)=ALL
OLOAD(NEU)=ALL
ELFO(NEU)=ALL
STRE(NEU)=ALL
STRN(NEU)=ALL
```

Expected behavior:

- `NEU` gets the requested vectors
- geometry snapshot is still written when `NEU` is active
- `F06` remains mainly diagnostic/summary text

Again, avoid bare `= ALL` if you do not want detailed `F06` tables.

### 3. `OP2` and `NEU`, but not detailed result tables in `F06`

Use:

```text
PARAM,OUTMODE,SMART
PARAM,PRTOP2,Y
PARAM,PRTF06,N
PARAM,PRTNEU,Y
```

Case Control examples:

```text
DISP(PLOT,NEU)=ALL
SPCF(PLOT,NEU)=ALL
OLOAD(PLOT,NEU)=ALL
ELFO(PLOT,NEU)=ALL
STRE(PLOT,NEU)=ALL
STRN(PLOT,NEU)=ALL
```

Expected behavior:

- `OP2` and `NEU` both receive output
- `F06` remains mainly summary/diagnostic

### 4. All three: detailed `F06` plus `OP2` plus `NEU`

Use either explicit descriptors:

```text
DISP(PRINT,PLOT,NEU)=ALL
SPCF(PRINT,PLOT,NEU)=ALL
OLOAD(PRINT,PLOT,NEU)=ALL
ELFO(PRINT,PLOT,NEU)=ALL
STRE(PRINT,PLOT,NEU)=ALL
STRN(PRINT,PLOT,NEU)=ALL
```

or broad flags:

```text
PARAM,PRTALL,Y
```

Expected behavior:

- detailed result tables in `F06`
- binary result payload in `OP2`
- neutral output in `NEU`

### 5. Only detailed `F06`

Use:

```text
PARAM,OUTMODE,LEGACY
PARAM,PRTF06,Y
PARAM,PRTOP2,N
PARAM,PRTNEU,N
```

or explicit Case Control:

```text
DISP(PRINT)=ALL
SPCF(PRINT)=ALL
OLOAD(PRINT)=ALL
ELFO(PRINT)=ALL
STRE(PRINT)=ALL
STRN(PRINT)=ALL
```

## Practical deck templates

### Template A: summary `F06` + `OP2` only

```text
SOL 101
CEND
TITLE = OP2 only test
SUBCASE 1
  DISP(PLOT)=ALL
  SPCF(PLOT)=ALL
  OLOAD(PLOT)=ALL
  STRE(PLOT)=ALL
  STRN(PLOT)=ALL
BEGIN BULK
PARAM,OUTMODE,SMART
PARAM,PRTOP2,Y
PARAM,PRTF06,N
PARAM,PRTNEU,N
ENDDATA
```

### Template B: summary `F06` + `NEU` only

```text
SOL 101
CEND
TITLE = NEU only test
SUBCASE 1
  DISP(NEU)=ALL
  SPCF(NEU)=ALL
  OLOAD(NEU)=ALL
  STRE(NEU)=ALL
  STRN(NEU)=ALL
BEGIN BULK
PARAM,OUTMODE,SMART
PARAM,PRTOP2,N
PARAM,PRTF06,N
PARAM,PRTNEU,Y
ENDDATA
```

### Template C: summary `F06` + `OP2` + `NEU`

```text
SOL 101
CEND
TITLE = OP2 and NEU without bulk F06 tables
SUBCASE 1
  DISP(PLOT,NEU)=ALL
  SPCF(PLOT,NEU)=ALL
  OLOAD(PLOT,NEU)=ALL
  STRE(PLOT,NEU)=ALL
  STRN(PLOT,NEU)=ALL
BEGIN BULK
PARAM,OUTMODE,SMART
PARAM,PRTOP2,Y
PARAM,PRTF06,N
PARAM,PRTNEU,Y
ENDDATA
```

### Template D: classic full output

```text
SOL 101
CEND
TITLE = full output
SUBCASE 1
  DISP(PRINT,PLOT,NEU)=ALL
  SPCF(PRINT,PLOT,NEU)=ALL
  OLOAD(PRINT,PLOT,NEU)=ALL
  STRE(PRINT,PLOT,NEU)=ALL
  STRN(PRINT,PLOT,NEU)=ALL
BEGIN BULK
PARAM,OUTMODE,SMART
PARAM,PRTOP2,Y
PARAM,PRTF06,Y
PARAM,PRTNEU,Y
ENDDATA
```

## Current limitations

This branch is in a transition state.

What is already true:

- `OUTMODE=SMART` exists
- `LINK9` no longer uses `PRTF06=Y` as a blanket “force all result detail” switch
- `NEU` routing is already split by family in `LINK9`

What is still not fully generalized:

- several Case Control family handlers still default bare `= ALL` to `PRINT+PLOT`
- so “`OP2 only`” or “`NEU only`” is reliable only when you use explicit descriptors
- `F06` is still always present as the main report/diagnostic file; this guide is about suppressing detailed result tables, not removing the file itself

## Source references

- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)
- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [LINK9.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/LINK9.f90)
- [CC_DISP.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_DISP.f90)
- [CC_STRE.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_STRE.f90)
- [CC_STRN.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_STRN.f90)
- [CC_ELFO.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_ELFO.f90)
- [CC_SPCF.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_SPCF.f90)
- [CC_OLOA.f90](D:/18a/MYSTRAN/Source/LK1/L1A-CC/CC_OLOA.f90)
