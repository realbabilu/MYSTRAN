# NEWSOLID continuum solid formulations

**Ary · Technical manuscript · 18 October 2026**

- [TETRA4 NEWSOLID: formulation and Nastran benchmarks](TETRA4_Formulation.md): nodal smoothing weight 0.9, mathematical operators, four figures and five existing verification cases.
- [Solids Formulation — read on GitHub](NEWSOLID_CONTINUUM_FORMULATIONS.md).
- [Alternative HTML version](NEWSOLID_CONTINUUM_FORMULATIONS.html): download this folder to view it locally with its figures. Equations are pre-rendered as SVG, with navigation and print styling. No network connection or equation-rendering script is required.
- Numerical tables and provenance (research data retained locally; see embedded tables).
- [Existing solid benchmark reports and input decks](../../solid_benchmark/README.md).

Eight solid topologies, with CPENTA6 ABE and EAS9 treated separately. Includes shape functions, explicit enhanced modes, condensation, smoothing, selective integration, pyramid quadrature, existing patch results, mesh figures and Nastran comparisons. All text is in English. The supplied revision-2 archive is reconciled with the actual MYSTRAN source snapshot; no new solver results are claimed.

Equations were generated using MathJax 3.2.2, distributed under Apache-2.0; see [license](assets/MathJax-LICENSE.txt).
