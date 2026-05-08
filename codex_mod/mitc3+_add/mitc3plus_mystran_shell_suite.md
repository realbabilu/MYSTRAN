# MITC3+ MYSTRAN Shell Suite

Generated from legacy `static-*mitc4*.bdf` shell decks by splitting each `CQUAD4` into two `CTRIA3`
and injecting `PARAM,TRIA3TYP,MITC3+`.

| Load Case | Theory | MYSTRAN MITC3+ | Error % |
|---|---:|---:|---:|
| `Static-22` | -3.09000000e-01 | -2.45034300e-01 | 20.701 |
| `Static-23` | 4.51970000e-04 | 4.56968000e-04 | 1.106 |
| `Static-24` | 9.40000000e-02 | 8.66592200e-02 | 7.809 |
| `Static-33 In-plane` | -5.42400000e-03 | -5.32377200e-03 | 1.848 |
| `Static-33 Out-of-plane` | -1.75400000e-03 | -1.87920400e-03 | 7.138 |
