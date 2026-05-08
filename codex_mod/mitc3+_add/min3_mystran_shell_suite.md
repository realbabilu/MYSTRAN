# MIN3 MYSTRAN Shell Suite

Generated from legacy `static-*mitc4*.bdf` shell decks by splitting each `CQUAD4` into two `CTRIA3`
and injecting `PARAM,TRIA3TYP,MIN3`.

| Load Case | Theory | MYSTRAN MIN3 | Error % |
|---|---:|---:|---:|
| `Static-22` | -3.09000000e-01 | -2.55462400e-01 | 17.326 |
| `Static-23` | 4.51970000e-04 | 4.60192900e-04 | 1.819 |
| `Static-24` | 9.40000000e-02 | 9.22110100e-02 | 1.903 |
| `Static-33 In-plane` | -5.42400000e-03 | -7.60020400e-03 | 40.122 |
| `Static-33 Out-of-plane` | -1.75400000e-03 | -2.20843900e-03 | 25.909 |
