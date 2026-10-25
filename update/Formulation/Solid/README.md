# NEWSOLID continuum solid formulations

[Continuum formulation overview](NEWSOLID_CONTINUUM_FORMULATIONS.md) covers all eight solid topologies, including interpolation, integration, enhanced modes, smoothing, recovery and existing Nastran benchmark comparisons. CTETRA10, CPENTA15, CHEXA20 and CPYRAM5 are covered in that overview, including the applicable selective-integration choices.

## Special element formulations

The following reports provide detailed mathematical formulations, figures and existing benchmark comparisons with Nastran:

| Element | Detailed report | Main formulation |
|---|---|---|
| TETRA4 | [TETRA4 NEWSOLID](TETRA4_Formulation.md) | Blended nodal strain smoothing, weight 0.9 |
| CPYRAM13 | [CPYRAM13 Bedrosian](CPYRAM13_Bedrosian_Formulation.md) | Rational quadratic pyramid interpolation |
| CHEXA8 | [CHEXA8 EAS](CHEXA8_EAS_Formulation.md) | Covariant EAS9 and consistent transverse-shear tying |
| CPENTA6 | [CPENTA6 ABE](CPENTA6_ABE_Formulation.md) | Edge-tied assumed shear and nine condensed enhanced strains |

## Benchmark collection

[Solid benchmark reports and input decks](../../solid_benchmark/README.md) provide the existing comparison results and model definitions. Equations and figures are included directly in the Markdown reports.
