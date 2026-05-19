! --- response_spectrum_mystran_add begin --- !
# FEMAP Output Profiles (PRINT / PLOT / NEU)

Tujuan: memudahkan pemilihan output MYSTRAN ala FEMAP.

## Mapping
- `PRINT` -> `F06` (text report)
- `PLOT` -> `OP2` (binary post file untuk FEMAP)
- `PRINT,PLOT` -> keduanya
- `NEU` -> bukan syntax Nastran standar; dibuat lewat post-step writer terpisah

## Profile 1: FAST (besar, cepat, SOL112-safe)
Gunakan `PLOT` saja untuk hasil detail RS yang stabil di SOL112.

```nastran
ECHO = NONE
SUBCASE 1
  METHOD = 20
  SPC = 10
  DLOAD = 30
  DISPLACEMENT(PLOT)=ALL
  VELOCITY(PLOT)=ALL
  ACCELERATION(PLOT)=ALL
  SPCFORCE(PLOT)=ALL
```

## Profile 2: CHECK (debug/validasi, SOL112-safe)
Gunakan `PRINT,PLOT` supaya ada trace lengkap di F06.

```nastran
ECHO = NONE
SUBCASE 1
  METHOD = 20
  SPC = 10
  DLOAD = 30
  DISPLACEMENT(PRINT,PLOT)=ALL
  VELOCITY(PRINT,PLOT)=ALL
  ACCELERATION(PRINT,PLOT)=ALL
  SPCFORCE(PRINT,PLOT)=ALL
```

## Profile 3: FEMAP (recommended workflow)
1. Solve MYSTRAN dengan `PLOT` sebagai sumber utama post.
2. Jika perlu file presentasi/import khusus FEMAP (`.neu`), generate via:
   - `scripts/femap_rs_neu_writer.py`

Contoh minimal (SOL112-safe):
```nastran
ECHO = NONE
SUBCASE 1
  METHOD = 20
  SPC = 10
  DLOAD = 30
  DISPLACEMENT(PLOT)=ALL
  VELOCITY(PLOT)=ALL
  ACCELERATION(PLOT)=ALL
  SPCFORCE(PLOT)=ALL
```

## P1 Profile Table (Exact Case Control)
| Profile | Tujuan | Exact Case Control | Output Utama |
|---|---|---|---|
| `FAST` | Runtime minimum untuk model besar | `ECHO=NONE`; `SUBCASE 1`; `METHOD=20`; `SPC=10`; `DLOAD=30`; request RS stabil: `DISPLACEMENT/VELOCITY/ACCELERATION/SPCFORCE (...PLOT)` | OP2 |
| `CHECK` | Verifikasi/debug hasil | `ECHO=NONE`; `SUBCASE 1`; `METHOD=20`; `SPC=10`; `DLOAD=30`; request RS stabil: `DISPLACEMENT/VELOCITY/ACCELERATION/SPCFORCE (...PRINT,PLOT)` | F06 + OP2 |
| `FEMAP` | Workflow FEMAP | Case Control sama dengan `FAST` (`PLOT`), lalu post-step `femap_rs_neu_writer.py` | OP2 + optional NEU |

Case Control template siap copy (P1):
```nastran
ECHO = NONE
SUBCASE 1
  METHOD = 20
  SPC = 10
  DLOAD = 30
  DISPLACEMENT(PLOT)=ALL
  VELOCITY(PLOT)=ALL
  ACCELERATION(PLOT)=ALL
  SPCFORCE(PLOT)=ALL
```

## Catatan Compatibility
- Untuk `SOL112` MYSTRAN saat fase ini, request `STRESS/STRAIN/FORCE` di profile RS dapat memicu error output-processing.
- Request tersebut dipindahkan ke fase lanjutan setelah jalur OFP/analysis code distabilkan.

## Catatan
- Untuk model besar, hindari `PRINT` penuh pada corner stress/strain kecuali benar-benar perlu.
- `NEU` diposisikan sebagai artefak post-processing, bukan solver-native output deck.

! --- response_spectrum_mystran_add end --- !
