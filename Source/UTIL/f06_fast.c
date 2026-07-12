#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static const char *DIGIT_TABLE =
    "00010203040506070809"
    "10111213141516171819"
    "20212223242526272829"
    "30313233343536373839"
    "40414243444546474849"
    "50515253545556575859"
    "60616263646566676869"
    "70717273747576777879"
    "80818283848586878889"
    "90919293949596979899";

static inline char *write2(char *dst, int val) {
    memcpy(dst, DIGIT_TABLE + val * 2, 2);
    return dst + 2;
}

static inline char *write6(char *dst, int val) {
    int hi = val / 10000;
    int lo = val % 10000;
    dst = write2(dst, hi);
    dst = write2(dst, lo / 100);
    dst = write2(dst, lo % 100);
    return dst;
}

static inline char *write_ndigits(char *dst, unsigned int val, int ndigits) {
    char buf[16];
    int i;

    for (i = ndigits - 1; i >= 0; --i) {
        buf[i] = (char)('0' + (val % 10u));
        val /= 10u;
    }

    memcpy(dst, buf, (size_t)ndigits);
    return dst + ndigits;
}

static inline int pow10i(int n) {
    static const int POWERS[] = {1, 10, 100, 1000, 10000, 100000, 1000000};
    if (n < 0 || n >= (int)(sizeof(POWERS) / sizeof(POWERS[0]))) {
        return 0;
    }
    return POWERS[n];
}

static inline char *fmt_i8_rj(char *dst, int val) {
    memset(dst, ' ', 8);

    if (val == 0) {
        dst[7] = '0';
        return dst + 8;
    }

    {
        int is_neg = (val < 0);
        unsigned int uval = is_neg ? (unsigned int)(-(val + 1)) + 1u : (unsigned int)val;
        int pos = 7;

        while (uval >= 100u) {
            unsigned int rem = uval % 100u;
            uval /= 100u;
            dst[pos - 1] = DIGIT_TABLE[rem * 2u];
            dst[pos]     = DIGIT_TABLE[rem * 2u + 1u];
            pos -= 2;
        }

        if (uval >= 10u) {
            dst[pos - 1] = DIGIT_TABLE[uval * 2u];
            dst[pos]     = DIGIT_TABLE[uval * 2u + 1u];
            pos -= 2;
        } else {
            dst[pos] = (char)('0' + uval);
            pos -= 1;
        }

        if (is_neg) {
            dst[pos] = '-';
        }
    }

    return dst + 8;
}

static char *fmt_e_width_prec(char *dst, double v, int width, int precision, int zero_special_f06) {
    int neg;
    int expv;
    int mant_int;
    int pad;
    int scale10;
    double mant;
    double scaled;

    if (zero_special_f06 && v == 0.0) {
        if (width == 14 && precision == 6) {
            memcpy(dst, "  0.0         ", 14);
            return dst + 14;
        }
        {
            int n = snprintf(dst, (size_t)width + 8u, "%*.*E", width, precision, v);
            return dst + n;
        }
    }

    if (v == 0.0 || isnan(v) || isinf(v)) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*E", width, precision, v);
        return dst + n;
    }

    neg = (v < 0.0);
    if (neg) {
        v = -v;
    }

    expv = (int)floor(log10(v));
    if (expv > 99 || expv < -99) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*E", width, precision, neg ? -v : v);
        return dst + n;
    }

    mant = v / pow(10.0, (double)expv);
    if (mant >= 10.0 - 1e-12) {
        mant /= 10.0;
        expv++;
    }
    if (mant < 1.0) {
        mant *= 10.0;
        expv--;
    }

    if (expv > 99 || expv < -99) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*E", width, precision, neg ? -v : v);
        return dst + n;
    }

    scale10 = pow10i(precision);
    if (scale10 <= 0) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*E", width, precision, neg ? -v : v);
        return dst + n;
    }

    scaled = mant * (double)scale10 + 0.5;
    mant_int = (int)scaled;
    if (mant_int >= 10 * scale10) {
        mant_int = scale10;
        expv++;
    }

    if (expv > 99 || expv < -99) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*E", width, precision, neg ? -v : v);
        return dst + n;
    }

    pad = width - (precision + 7);
    while (pad-- > 0) {
        *dst++ = ' ';
    }

    if (neg) {
        *dst++ = '-';
    } else {
        *dst++ = ' ';
    }

    *dst++ = (char)('0' + (mant_int / scale10));
    *dst++ = '.';
    if (precision == 6) {
        dst = write6(dst, mant_int % scale10);
    } else {
        dst = write_ndigits(dst, (unsigned int)(mant_int % scale10), precision);
    }
    *dst++ = 'E';
    if (expv >= 0) {
        *dst++ = '+';
    } else {
        *dst++ = '-';
        expv = -expv;
    }
    dst = write2(dst, expv);
    return dst;
}

static char *fmt_f8_2(char *dst, double v) {
    int neg;
    int pad;
    uint64_t scaled;
    uint64_t intpart;
    unsigned int frac;
    char ibuf[24];
    int ilen;

    if (isnan(v) || isinf(v)) {
        int n = snprintf(dst, 16u, "%8.2f", v);
        return dst + n;
    }

    neg = (v < 0.0);
    if (neg) {
        v = -v;
    }

    scaled = (uint64_t)(v * 100.0 + 0.5);
    intpart = scaled / 100u;
    frac = (unsigned int)(scaled % 100u);

    ilen = 0;
    do {
        ibuf[ilen++] = (char)('0' + (int)(intpart % 10u));
        intpart /= 10u;
    } while (intpart != 0u && ilen < (int)sizeof(ibuf));

    pad = 8 - (ilen + 1 + 2 + (neg ? 1 : 0));
    if (pad < 0) {
        int n = snprintf(dst, 16u, "%8.2f", neg ? -v : v);
        return dst + n;
    }

    while (pad-- > 0) {
        *dst++ = ' ';
    }
    if (neg) {
        *dst++ = '-';
    }
    while (ilen-- > 0) {
        *dst++ = ibuf[ilen];
    }
    *dst++ = '.';
    dst = write2(dst, (int)frac);
    return dst;
}

static char *fmt_f_width_prec(char *dst, double v, int width, int precision) {
    int neg;
    int pad;
    int scale10;
    uint64_t scaled;
    uint64_t intpart;
    unsigned int frac;
    char ibuf[32];
    int ilen;

    if (isnan(v) || isinf(v)) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*f", width, precision, v);
        return dst + n;
    }

    neg = (v < 0.0);
    if (neg) {
        v = -v;
    }

    scale10 = pow10i(precision);
    if (scale10 <= 0) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*f", width, precision, neg ? -v : v);
        return dst + n;
    }

    scaled = (uint64_t)(v * (double)scale10 + 0.5);
    intpart = scaled / (uint64_t)scale10;
    frac = (unsigned int)(scaled % (uint64_t)scale10);

    ilen = 0;
    do {
        ibuf[ilen++] = (char)('0' + (int)(intpart % 10u));
        intpart /= 10u;
    } while (intpart != 0u && ilen < (int)sizeof(ibuf));

    pad = width - (ilen + 1 + precision + (neg ? 1 : 0));
    if (pad < 0) {
        int n = snprintf(dst, (size_t)width + 8u, "%*.*f", width, precision, neg ? -v : v);
        return dst + n;
    }

    while (pad-- > 0) {
        *dst++ = ' ';
    }
    if (neg) {
        *dst++ = '-';
    }
    while (ilen-- > 0) {
        *dst++ = ibuf[ilen];
    }
    *dst++ = '.';
    dst = write_ndigits(dst, frac, precision);
    return dst;
}

void mystran_fmt_e14_6_f06(double v, char out[14]) {
    (void)fmt_e_width_prec(out, v, 14, 6, 1);
}

void mystran_fmt_e17_6(double v, char out[17]) {
    (void)fmt_e_width_prec(out, v, 17, 6, 0);
}

void mystran_fmt_es13_5(double v, char out[13]) {
    (void)fmt_e_width_prec(out, v, 13, 5, 0);
}

void mystran_fmt_es14_5(double v, char out[14]) {
    (void)fmt_e_width_prec(out, v, 14, 5, 0);
}

void mystran_fmt_es11_3(double v, char out[11]) {
    (void)fmt_e_width_prec(out, v, 11, 3, 0);
}

void mystran_fmt_es10_2(double v, char out[10]) {
    (void)fmt_e_width_prec(out, v, 10, 2, 0);
}

void mystran_fmt_f8_2(double v, char out[8]) {
    (void)fmt_f8_2(out, v);
}

void mystran_fmt_f9_3(double v, char out[9]) {
    (void)fmt_f_width_prec(out, v, 9, 3);
}

void mystran_fmt_e9_1(double v, char out[9]) {
    (void)fmt_e_width_prec(out, v, 9, 1, 0);
}

void mystran_fmt_i8_rj(int val, char out[8]) {
    (void)fmt_i8_rj(out, val);
}

void mystran_build_grid_f06_line(int gid, int coord, const double vals[6], char out[108]) {
    int j;
    char *p = out;

    memcpy(p, "      ", 6);
    p += 6;
    *p++ = ' ';
    p = fmt_i8_rj(p, gid);
    *p++ = ' ';
    p = fmt_i8_rj(p, coord);
    for (j = 0; j < 6; ++j) {
        p = fmt_e_width_prec(p, vals[j], 14, 6, 1);
    }
}

void mystran_build_quad_1403_line(int eid, const double vals[10], char out[145]) {
    int j;
    char *p = out;

    *p++ = ' ';
    p = fmt_i8_rj(p, eid);
    memcpy(p, "  CENTER     ", 13);
    p += 13;
    p = fmt_e_width_prec(p, vals[0], 11, 3, 0);
    for (j = 1; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_f8_2(p, vals[4]);
    for (j = 5; j < 10; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
}

void mystran_build_quad_1404_line(const double vals[8], char out[119]) {
    int j;
    char *p = out;

    memset(p, ' ', 22);
    p += 22;
    p = fmt_e_width_prec(p, vals[0], 11, 3, 0);
    for (j = 1; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_f8_2(p, vals[4]);
    for (j = 5; j < 8; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
}

void mystran_build_quad_1405_line(int gid, const double vals[10], double poly_err, int poly_idx, char out[157]) {
    int j;
    char *p = out;

    memset(p, ' ', 11);
    p += 11;
    memcpy(p, "GRD", 3);
    p += 3;
    p = fmt_i8_rj(p, gid);
    p = fmt_e_width_prec(p, vals[0], 11, 3, 0);
    for (j = 1; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_f8_2(p, vals[4]);
    for (j = 5; j < 10; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_e_width_prec(p, poly_err, 9, 1, 0);
    *p++ = '(';
    *p++ = (char)('0' + poly_idx);
    *p++ = ')';
}

void mystran_build_quad_1406_line(int gid, const double vals[10], double poly_err, char out[154]) {
    int j;
    char *p = out;

    memset(p, ' ', 11);
    p += 11;
    memcpy(p, "GRD", 3);
    p += 3;
    p = fmt_i8_rj(p, gid);
    p = fmt_e_width_prec(p, vals[0], 11, 3, 0);
    for (j = 1; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_f8_2(p, vals[4]);
    for (j = 5; j < 10; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_e_width_prec(p, poly_err, 9, 1, 0);
}

void mystran_build_tria_1703_line(int eid, const double vals[10], char out[149]) {
    int j;
    char *p = out;

    *p++ = ' ';
    p = fmt_i8_rj(p, eid);
    memcpy(p, "    Anywhere  ", 14);
    p += 14;
    for (j = 0; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_e_width_prec(p, vals[4], 9, 3, 0);
    for (j = 5; j < 10; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
}

void mystran_build_tria_1704_line(const double vals[8], char out[149]) {
    int j;
    char *p = out;

    memset(p, ' ', 13);
    p += 13;
    memcpy(p, "in elem", 7);
    p += 7;
    memset(p, ' ', 3);
    p += 3;
    for (j = 0; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_e_width_prec(p, vals[4], 9, 3, 0);
    for (j = 5; j < 8; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
}

void mystran_build_tria_1706_line(int gid, const double vals[10], char out[139]) {
    int j;
    char *p = out;

    *p++ = ' ';
    p = fmt_i8_rj(p, gid);
    memset(p, ' ', 4);
    p += 4;
    for (j = 0; j <= 3; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
    p = fmt_e_width_prec(p, vals[4], 9, 3, 0);
    for (j = 5; j < 10; ++j) {
        p = fmt_e_width_prec(p, vals[j], 13, 5, 0);
    }
}
