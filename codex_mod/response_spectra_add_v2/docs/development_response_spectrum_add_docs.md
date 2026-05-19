! --- response_spectrum_mystran_add begin --- !
# Development Notes

## Architecture
- Keep RS neutral generation outside MYSTRAN core.
- Use existing MYSTRAN outputs as source:
  - SOL103 NEU/F06 for modes and frequencies
  - SOL112 NEU for RSX/RSY vectors
- Compose final FEMAP NEU with geometry + output sets.

## Why External
- Faster iteration for FEMAP-neutral formatting.
- Avoid touching core solver while RS UI/export still evolving.
- Reduce risk around proprietary neutral interpretation concerns.

## Combo Strategy
- Base directional sets: RSX, RSY
- Sign combos: +/? variants
- Civil presets: 1.0X+0.3Y and 0.3X+1.0Y (+ sign variants)
- All computed by per-node linear superposition in writer.
- Writer now supports CLI-defined combos:
  - `--combo "NAME:a:b"` => `a*RSX + b*RSY`
  - `--no-default-combos` to disable built-in combo presets.

## Current Limitation
- This layer does not replace solver-native CQC3 directional combination.
- It is a deterministic linear post-combination layer.

! --- response_spectrum_mystran_add end --- !
