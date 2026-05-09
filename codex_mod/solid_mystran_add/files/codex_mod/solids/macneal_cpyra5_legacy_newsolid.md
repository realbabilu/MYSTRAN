# MacNeal Static-37 CPYRA5 Legacy vs NEWSOLID

Same CPYRA5 pyramid mesh; legacy means standard PYRA5, NEWSOLID means CPYRA5_EAS54 stiffness.

| mesh | case | reference | legacy avg | legacy/ref | NEWSOLID avg | NEWSOLID/ref | NEWSOLID/legacy |
|---|---|---:|---:|---:|---:|---:|---:|
| `4x2x1 blocks -> 6 pyra/block` | `inplane_z` | `-5.424000000e-03` | `-1.207589667e-04` | `0.022` | `-1.837215667e-04` | `0.034` | `1.521` |
| `4x2x1 blocks -> 6 pyra/block` | `outplane_y` | `-1.754000000e-03` | `-9.680725000e-05` | `0.055` | `-1.396055000e-04` | `0.080` | `1.442` |
| `8x2x1 blocks -> 6 pyra/block` | `inplane_z` | `-5.424000000e-03` | `-4.393531000e-04` | `0.081` | `-6.739854333e-04` | `0.124` | `1.534` |
| `8x2x1 blocks -> 6 pyra/block` | `outplane_y` | `-1.754000000e-03` | `-2.868131333e-04` | `0.164` | `-4.110636000e-04` | `0.234` | `1.433` |
| `12x2x1 blocks -> 6 pyra/block` | `inplane_z` | `-5.424000000e-03` | `-8.694420333e-04` | `0.160` | `-1.349661000e-03` | `0.249` | `1.552` |
| `12x2x1 blocks -> 6 pyra/block` | `outplane_y` | `-1.754000000e-03` | `-4.693927000e-04` | `0.268` | `-6.753357333e-04` | `0.385` | `1.439` |
