# External MUMPS folder placeholder

This `mumps` folder is intentionally left without the real `MUMPS` source or
library payload.

Reason:

- `MUMPS` is an external dependency
- users may or may not want to enable it
- this patch bundle is designed so the core source changes can be uploaded
  without vendoring the external dependency itself

## Expected user action

If `MUMPS` support is desired, the user should manually download/build
sequential non-MPI `MUMPS` and place the resulting source/build tree under the
repository `mumps/` folder.

Recommended bootstrap helper:

- [scivision/mumps-superbuild](https://github.com/scivision/mumps-superbuild)

See also:

- `..\..\mumps_nonmpi_static_library_handoff.md`

## Important

Do not expect this placeholder folder alone to enable `MUMPS`.
