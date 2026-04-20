# CBEAM Modifications Overview (MYSTRAN 17)

This document summarizes the CBEAM work added in this branch and how it flows through MYSTRAN.

## Scope

Implemented and validated:

- 3D, 2-node, 12-DOF `CBEAM` stiffness path in `BEAM.f90`
- DSB/Timoshenko-style bending behavior in both bending planes via internal shear influence factors
- axial + torsion + two bending planes in one local 12x12 element matrix
- thermal load/recovery reuse through existing 1D thermal pipeline
- end release support via existing `DOFPIN -> PINFLG` mechanism
- `PLOAD1` equivalent nodal loads for:
  - local y (`FYE/FY/Y`)
  - local z (`FZE/FZ/Z`)
  - local x-force (`FXE`)
  - local x-moment (`MXE`)
  - local y-moment (`MYE`)
  - local z-moment (`MZE`)
  - full-span, partial-span, and concentrated-in-element (`X1=X2`)
- engineering-force + 1D stress/strain recovery path no longer hard-fails for BEAM

## Main Source Files Touched

- `Source/EMG/EMG3/BEAM.f90`
- `Source/EMG/EMG3/BAR1.f90` (same `PLOAD1` handling pattern for BAR)
- `Source/LK1/L1D/PRESSURE_DATA_PROC.f90`
- `Source/EMG/EMG1/ELMDAT2.f90`
- `Source/Modules/SCONTR.f90`
- `Source/Modules/MODEL_STUF.f90`
- `Source/LK1/L1A-BD/BD_PLOAD2.f90`
- `Source/LK9/L92/OFP3_ELFE_1D.f90`
- `Source/LK9/L92/ONE_D_STRESS_OUTPUTS.f90`
- `Source/LK9/L92/ONE_D_STRAIN_OUTPUTS.f90`

## CBEAM Element Formulation Flow (`BEAM.f90`)

```mermaid
flowchart TD
  A["ELMDAT1 provides A,I1,I2,I12,J,K1,K2,E,G,ALPHA,TREF,L"] --> B["Compute phi1, phi2 internally"]
  B --> C["Build KE blocks:
  axial (UX),
  torsion (RX),
  bending plane-1 (UY,RZ),
  bending plane-2 (UZ,RY)"]
  C --> D["Symmetrize KE"]
  D --> E["Apply end releases via DOFPIN/PINFLG"]
  E --> F["Build thermal quantities TPRIME, ABAR"]
  F --> G["Compute PTE and STE using existing 1D convention"]
  G --> H["If nonlinear step: build KED from current element forces"]
```

## `PLOAD1` Pipeline Flow (Data To Solve)

```mermaid
flowchart TD
  A["Bulk Data: PLOAD1"] --> B["BD_PLOAD2.f90
  accepts beam local types + FR scale"]
  B --> C["PRESSURE_DATA_PROC.f90
  validates X1/X2 and stores per component:
  [P1,P2,X1,X2] x 4"]
  C --> D["PPNT/PDATA/PTYPE written"]
  D --> E["ELMDAT2.f90 copies 16-slot PRESS for beam/bar"]
  E --> F["BEAM.f90 / BAR1.f90
  convert PRESS to PPE"]
  F --> G["EMG assembly:
  PPE -> EPTL -> SYS_LOAD -> PG"]
  G --> H["Link2 reductions -> Link3 solve"]
  H --> I["Link9 recovery + output"]
```

## `PLOAD1` Equivalent Nodal Conversion Logic (`BEAM.f90`)

```mermaid
flowchart TD
  A["Read one component block [P1,P2,X1,X2]"] --> B{"X1 < 0 ?"}
  B -- Yes --> C["Component absent: skip"]
  B -- No --> D{"X2 == X1 ?"}
  D -- Yes --> E["Concentrated in element:
  evaluate shape functions at x=X1"]
  D -- No --> F["Distributed on [X1,X2]:
  3-point Gauss integrate shape functions"]
  E --> G["Accumulate into PPE DOFs"]
  F --> G
  G --> H["Repeat for y, z, x-force, x-moment components"]
```

## Recovery/Output Flow (Static `SOL 101`)

```mermaid
flowchart TD
  A["Solved displacements UG"] --> B["Element recovery builds PEL/SE/STE"]
  B --> C["OFP3_ELFE_1D.f90 maps BEAM engineering forces"]
  B --> D["ONE_D_STRESS_OUTPUTS.f90 includes BEAM"]
  B --> E["ONE_D_STRAIN_OUTPUTS.f90 includes BEAM"]
  C --> F["F06 engineering-force table"]
  D --> G["F06 stress table"]
  E --> H["F06 strain table"]
```

## Internal `phi` Treatment

`phi` is not a new input field in v1. It is derived internally per element:

- `phi1 = 12*E*I1 / (K2*G*A*L^2)` for plane-1 block (`UY/RZ`)
- `phi2 = 12*E*I2 / (K1*G*A*L^2)` for plane-2 block (`UZ/RY`)

This keeps existing `PBEAM` inputs unchanged while enabling Timoshenko/DSB-like shear flexibility behavior.

## Current `PLOAD1` Behavior (Beam/Bar)

Supported now:

- full-span uniform load
- partial-span distributed load
- concentrated-in-element load (`X1=X2`)
- linear variation along loaded segment (`P1` to `P2`)

Current constraints:

- `SCALE=FR` only
- one `PLOAD1` entry per `(subcase, element, component)` in current implementation
- no global-direction `FX/FY/FZ/MX/MY/MZ` beam interpretation in this path

## Validation Decks

- [cbeam_validation_cantilever_point_load.md](E:/mystran17/mystran/dev_docs/cbeam_validation_cantilever_point_load.md)
- [cbeam_validation_axial_extension.md](E:/mystran17/mystran/dev_docs/cbeam_validation_axial_extension.md)
- [cbeam_validation_torsion.md](E:/mystran17/mystran/dev_docs/cbeam_validation_torsion.md)
- [cbeam_validation_portal_frame_2d.md](E:/mystran17/mystran/dev_docs/cbeam_validation_portal_frame_2d.md)
- [cbeam_validation_pload1.md](E:/mystran17/mystran/dev_docs/cbeam_validation_pload1.md)
- [cbeam_validation_pload1_axial_torsion.md](E:/mystran17/mystran/dev_docs/cbeam_validation_pload1_axial_torsion.md)
- [cbeam_validation_pload1_partial_concentrated.md](E:/mystran17/mystran/dev_docs/cbeam_validation_pload1_partial_concentrated.md)
- [cbeam_validation_pload1_mye_mze.md](E:/mystran17/mystran/dev_docs/cbeam_validation_pload1_mye_mze.md)

## Notes For GitHub Push

- This file is intended as the top-level summary for CBEAM changes.
- Mermaid diagrams render directly on GitHub in markdown view.
- If you want, we can next split this into `LINK1`, `EMG`, and `LINK9` sub-pages and link them as a mini index.
