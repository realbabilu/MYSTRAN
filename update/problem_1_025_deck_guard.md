# Problem 1-025 Deck Guard

Date: July 14, 2026

## Correct interpretation

- `D:\18a\femap_RSA\problem_1_025_rsa_srss_new.dat`
  is the active MYSTRAN response-spectrum deck for `SAP2000 Problem 1-025`.
- The correct modal basis check for that deck must come from the same structural model with the RSA application cards removed.

## Wrong comparison deck for RSA

- `D:\18a\femap_RSA\problem_1_025_mystran_dense_modal_axial_surrogate.dat`
  is an older axial-surrogate modal package.
- It is useful only as a separate surrogate study.
- It must not be used as the modal reference for:
  - `problem_1_025_rsa_srss_new.dat`
  - `problem_1_025_rsa.dat`

Reason:

- the surrogate deck is not the same structural model or property set as the current RSA deck
- therefore mismatched eigenvalues between:
  - `problem_1_025_mystran_dense_modal_axial_surrogate.*`
  - `problem_1_025_rsa_srss_new.*`
  do not indicate an RSA solver bug

## Correct local MYSTRAN check created in this pass

- Modal-only deck derived directly from the active RSA deck:
  - `D:\18a\femap_RSA\problem_1_025_modal_from_rsa_dense.dat`

Expected result:

- `problem_1_025_modal_from_rsa_dense.*` and `problem_1_025_rsa_srss_new.*`
  should share the same eigen basis before RSA combination.

Observed basis:

- mode 1: `3.662793 Hz`
- mode 2: `6.342935 Hz`

## Operational rule

- For MYSTRAN-only debugging of `1-025`, compare RSA against:
  - `problem_1_025_modal_from_rsa_dense.*`
- Do not compare RSA against:
  - `problem_1_025_mystran_dense_modal_axial_surrogate.*`
