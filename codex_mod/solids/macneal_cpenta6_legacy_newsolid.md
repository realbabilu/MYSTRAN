# MacNeal Static-37 CPENTA6 Legacy vs NEWSOLID

Same CPENTA6 wedge mesh; legacy means standard CPENTA6, NEWSOLID means CPENTA6_EAS9.

| mesh | case | reference | legacy avg | legacy/ref | NEWSOLID avg | NEWSOLID/ref | NEWSOLID/legacy |
|---|---|---:|---:|---:|---:|---:|---:|
| `4x2x1 blocks -> 2 penta/block` | `inplane_z` | `-5.424000000e-03` | `-1.503850833e-04` | `0.028` | `-2.855652000e-04` | `0.053` | `1.899` |
| `4x2x1 blocks -> 2 penta/block` | `outplane_y` | `-1.754000000e-03` | `-7.551971500e-05` | `0.043` | `-1.237743167e-04` | `0.071` | `1.639` |
| `8x2x1 blocks -> 2 penta/block` | `inplane_z` | `-5.424000000e-03` | `-5.712354667e-04` | `0.105` | `-1.010580167e-03` | `0.186` | `1.769` |
| `8x2x1 blocks -> 2 penta/block` | `outplane_y` | `-1.754000000e-03` | `-2.503923833e-04` | `0.143` | `-4.074410167e-04` | `0.232` | `1.627` |
| `12x2x1 blocks -> 2 penta/block` | `inplane_z` | `-5.424000000e-03` | `-1.151547167e-03` | `0.212` | `-1.908320833e-03` | `0.352` | `1.657` |
| `12x2x1 blocks -> 2 penta/block` | `outplane_y` | `-1.754000000e-03` | `-4.472169500e-04` | `0.255` | `-7.073843833e-04` | `0.403` | `1.582` |
