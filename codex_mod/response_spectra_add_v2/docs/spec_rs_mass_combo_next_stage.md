! --- response_spectrum_mystran_add begin --- !
# Next Stage Spec: RS Mass Source + Combination (P1 Extension)

## Goal
Menambahkan kemampuan dasar ala workflow sipil:
1. Load-to-Mass (selain self weight)
2. Self weight sebagai sumber massa dinamik
3. Kombinasi hasil: Static+Static, Static+RS dengan envelope max/min

## Scope
- Fokus: implementasi internal + output konsisten untuk post-processing.
- Non-scope tahap ini:
  - Nonlinear transient
  - CQC3 multi-direction advanced
  - Native FEMAP integrated generator di core solver

## 1) Mass Source Control
### Konsep keyword internal
- `RSMASS, SID, MATMASS, YES/NO, SWMASS, YES/NO, G, <value>, DIR, X/Y/Z`
  - `MATMASS`: gunakan density material sebagai sumber massa.
  - `SWMASS`: aktifkan self-weight sebagai sumber massa ekuivalen.
  - `G`: percepatan gravitasi referensi konversi beban->massa.
  - `DIR`: arah gravitasi global.

### Load-to-Mass tambahan
- `RSMASSL, SID, LOADSET, LID, SCALE, s`
- Konversi: `m_eq = SCALE * W/g`

### Rule anti double-count
- Jika `MATMASS=YES` dan `SWMASS=YES`, solver wajib:
  1. mencatat sumber kontribusi,
  2. mencegah input yang menduplikasi massa yang sama,
  3. menulis warning eksplisit di log.

## 2) Self Weight as Dynamic Mass
Jika `SWMASS=YES`:
1. ambil vektor gaya gravitasi aktif,
2. konversi ke kontribusi massa dengan `g` referensi,
3. akumulasi ke matriks massa dinamik.

Kriteria:
- massa total terlapor berubah sesuai aktivasi flag.
- mode shape/frequency berubah konsisten saat SWMASS on/off.

## 3) Combination Engine
### Konsep keyword internal
- `RSCOMB, CID, TYPE, LINEAR/ENVELOPE`
- `RSCOMBI, CID, TERM, CASEID, FACTOR`
- `RSCOMBRS, CID, RSCASE, METHOD, SRSS/CQC, SIGN, PM/MP/PP/MM`

### Capability minimum
1. `Static + Static` (linear superposition)
2. `Static + RS` (linear) + envelope `MAX/MIN`

### Output envelope minimum
- Displacement: `T1,T2,T3,R1,R2,R3`
- Reactions / constraint forces
- Element force/stress resultants

## 4) Data Flow
1. Build mass from configured sources (`MATMASS`, `SWMASS`, `RSMASSL`).
2. Solve modal base.
3. Build RS directional results.
4. Build RS combos (sign + 100/30 as requested set).
5. Combine with static cases.
6. Emit `MAX/MIN` envelopes.

## 5) Validation Plan
### Bench Set
- RS-01 (single direction baseline)
- RS-04 (X/Y directional + combo)
- One static gravity + lateral static case for mixed combination

### Acceptance checks
1. `MATMASS only` run valid.
2. `MATMASS + SWMASS` run valid.
3. `LOADSET->MASS` contribution terlapor dan konsisten.
4. `Static+Static` nilai cocok hitung manual.
5. `Static+RS` punya output `MAX/MIN` untuk displacement + force.

## 6) Suggested Implementation Order
1. Parser/keyword internal mass-source control.
2. Mass assembly hooks + reporting.
3. Static+Static combiner.
4. Static+RS combiner + envelope.
5. Validation scripts + sample artifacts.

## 7) Risks
1. Double-count massa dari source campuran.
2. Konvensi tanda saat envelope lintas direction/combination.
3. Konsistensi unit pada konversi `W/g`.

## 8) Deliverables
- Spec ini (approved)
- Contoh input P1 extension
- Validation report `Finished solution / Failed validation`
- Updated user guide section for mass-source + combination

! --- response_spectrum_mystran_add end --- !
