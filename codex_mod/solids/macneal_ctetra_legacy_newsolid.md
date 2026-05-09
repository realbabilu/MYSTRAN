# MacNeal Static-37 CTETRA Legacy vs NEWSOLID

Same tetra meshes; CTETRA4 NEWSOLID is smooth alpha 0.9, CTETRA10 NEWSOLID is the baseline quadratic guard.

| type | mesh | case | reference | legacy avg | legacy/ref | NEWSOLID avg | NEWSOLID/ref | NEWSOLID/legacy |
|---|---|---|---:|---:|---:|---:|---:|---:|
| `CTETRA4` | `4x2x1 blocks -> 6 tet/block` | `inplane_z` | `-5.424000000e-03` | `-6.872046000e-05` | `0.013` | `-4.832198333e-04` | `0.089` | `7.032` |
| `CTETRA4` | `4x2x1 blocks -> 6 tet/block` | `outplane_y` | `-1.754000000e-03` | `-4.638082000e-05` | `0.026` | `-2.942085000e-04` | `0.168` | `6.343` |
| `CTETRA4` | `8x2x1 blocks -> 6 tet/block` | `inplane_z` | `-5.424000000e-03` | `-2.437181000e-04` | `0.045` | `-1.795422667e-03` | `0.331` | `7.367` |
| `CTETRA4` | `8x2x1 blocks -> 6 tet/block` | `outplane_y` | `-1.754000000e-03` | `-1.444362667e-04` | `0.082` | `-8.035865333e-04` | `0.458` | `5.564` |
| `CTETRA4` | `12x2x1 blocks -> 6 tet/block` | `inplane_z` | `-5.424000000e-03` | `-4.794252000e-04` | `0.088` | `-3.407183333e-03` | `0.628` | `7.107` |
| `CTETRA4` | `12x2x1 blocks -> 6 tet/block` | `outplane_y` | `-1.754000000e-03` | `-2.549601000e-04` | `0.145` | `-1.311669333e-03` | `0.748` | `5.145` |
| `CTETRA10` | `4x2x1 blocks -> 6 tet/block` | `inplane_z` | `-5.424000000e-03` | `-3.083639000e-03` | `0.569` | `-3.083639000e-03` | `0.569` | `1.000` |
| `CTETRA10` | `4x2x1 blocks -> 6 tet/block` | `outplane_y` | `-1.754000000e-03` | `-1.082875667e-03` | `0.617` | `-1.082875667e-03` | `0.617` | `1.000` |
| `CTETRA10` | `8x2x1 blocks -> 6 tet/block` | `inplane_z` | `-5.424000000e-03` | `-3.862788000e-03` | `0.712` | `-3.862788000e-03` | `0.712` | `1.000` |
| `CTETRA10` | `8x2x1 blocks -> 6 tet/block` | `outplane_y` | `-1.754000000e-03` | `-1.285705000e-03` | `0.733` | `-1.285705000e-03` | `0.733` | `1.000` |
