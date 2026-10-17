#include <locale.h>

/* BDF real fields use a decimal point regardless of the Windows locale.
 * Reassert immediately before conversion: opening internal unformatted
 * files can restore the process locale in the MinGW Fortran runtime.
 */
void mystran_set_numeric_locale_(void) {
    setlocale(LC_NUMERIC, "C");
}
