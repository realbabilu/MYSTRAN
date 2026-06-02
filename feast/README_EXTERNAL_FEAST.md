# External FEAST folder placeholder

This `feast` folder is intentionally left without the FEAST source payload.

Reason:
- FEAST v4.0 is an external dependency
- users may or may not want to enable it
- the MYSTRAN source patch is designed so FEAST is optional

## Expected user action

If FEAST support is desired, the user should manually download and unpack FEAST
v4.0 into the repository's `feast` folder after uploading this patch.

Typical intent:

- repository root
  - `feast/`
    - FEAST source tree

The CMake and solver-side patch are written so:

- builds still work without FEAST when FEAST is not enabled
- FEAST paths become active only when the external FEAST library is provided and
  the corresponding CMake option is turned on

## Important

Do not expect this placeholder folder alone to enable FEAST.
