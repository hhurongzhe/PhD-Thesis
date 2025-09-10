#pragma once
#ifndef APWP_HPP
#define APWP_HPP

#include <iostream>
#include <cmath>
#include <tuple>
#include <complex>

namespace apwp
{
typedef std::complex<double> cdouble;

inline int32_t compute_key(int m1, int m2, int m3, int m4) { return ((m1 + 1) / 2) | (((m2 + 1) / 2) << 1) | (((m3 + 1) / 2) << 2) | (((m4 + 1) / 2) << 3); }

inline std::tuple<int, int, int, int> decode_key(int32_t key)
{
    int m1 = 2 * (key & 1) - 1;
    int m2 = 2 * ((key >> 1) & 1) - 1;
    int m3 = 2 * ((key >> 2) & 1) - 1;
    int m4 = 2 * ((key >> 3) & 1) - 1;
    return {m1, m2, m3, m4};
}

inline double Power(double x, int n) { return std::pow(x, n); }

inline cdouble apwd(const double &f1, const double &f2, const double &f3, const double &f4, const double &f5, const double &f6, const double &pfx, const double &pfy, const double &pfz, const double &pix, const double &piy, const double &piz, const int &m1, const int &m2, const int &m3, const int &m4)
{
    double rep = 0.0;
    double imp = 0.0;
    int32_t key = compute_key(m1, m2, m3, m4);
    switch (key)
    {
        case 0:
            rep = f1 + f2 + f3 * Power(pfz, 2) + (f4 * Power(pfz, 2)) / 4. + f6 * Power(pfy, 2) * Power(pix, 2) - 2 * f6 * pfx * pfy * pix * piy + f6 * Power(pfx, 2) * Power(piy, 2) - 2 * f3 * pfz * piz + (f4 * pfz * piz) / 2. + f3 * Power(piz, 2) + (f4 * Power(piz, 2)) / 4.;
            imp = -(f5 * pfy * pix) + f5 * pfx * piy;
            return cdouble(rep, imp);
        case 1:
            rep = -(f3 * pfx * pfz) - (f4 * pfx * pfz) / 4. + f3 * pfz * pix - (f4 * pfz * pix) / 4. - (f5 * pfz * pix) / 2. - f6 * pfy * pfz * pix * piy + f6 * pfx * pfz * Power(piy, 2) + f3 * pfx * piz - (f4 * pfx * piz) / 4. + (f5 * pfx * piz) / 2. - f3 * pix * piz - (f4 * pix * piz) / 4. +
                  f6 * Power(pfy, 2) * pix * piz - f6 * pfx * pfy * piy * piz;
            imp = f3 * pfy * pfz + (f4 * pfy * pfz) / 4. - f6 * pfy * pfz * Power(pix, 2) - f3 * pfz * piy + (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. + f6 * pfx * pfz * pix * piy - f3 * pfy * piz + (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. + f6 * pfx * pfy * pix * piz + f3 * piy * piz + (f4 * piy * piz) / 4. -
                  f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 2:
            rep = -(f3 * pfx * pfz) - (f4 * pfx * pfz) / 4. + f3 * pfz * pix - (f4 * pfz * pix) / 4. - (f5 * pfz * pix) / 2. - f6 * pfy * pfz * pix * piy + f6 * pfx * pfz * Power(piy, 2) + f3 * pfx * piz - (f4 * pfx * piz) / 4. + (f5 * pfx * piz) / 2. - f3 * pix * piz - (f4 * pix * piz) / 4. +
                  f6 * Power(pfy, 2) * pix * piz - f6 * pfx * pfy * piy * piz;
            imp = f3 * pfy * pfz + (f4 * pfy * pfz) / 4. - f6 * pfy * pfz * Power(pix, 2) - f3 * pfz * piy + (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. + f6 * pfx * pfz * pix * piy - f3 * pfy * piz + (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. + f6 * pfx * pfy * pix * piz + f3 * piy * piz + (f4 * piy * piz) / 4. -
                  f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 3:
            rep = f3 * Power(pfx, 2) + (f4 * Power(pfx, 2)) / 4. - f3 * Power(pfy, 2) - (f4 * Power(pfy, 2)) / 4. - 2 * f3 * pfx * pix + (f4 * pfx * pix) / 2. + f3 * Power(pix, 2) + (f4 * Power(pix, 2)) / 4. - f6 * Power(pfz, 2) * Power(pix, 2) + 2 * f3 * pfy * piy - (f4 * pfy * piy) / 2. - f3 * Power(piy, 2) -
                  (f4 * Power(piy, 2)) / 4. + f6 * Power(pfz, 2) * Power(piy, 2) + 2 * f6 * pfx * pfz * pix * piz - 2 * f6 * pfy * pfz * piy * piz - f6 * Power(pfx, 2) * Power(piz, 2) + f6 * Power(pfy, 2) * Power(piz, 2);
            imp = -2 * f3 * pfx * pfy - (f4 * pfx * pfy) / 2. + 2 * f3 * pfy * pix - (f4 * pfy * pix) / 2. + 2 * f3 * pfx * piy - (f4 * pfx * piy) / 2. - 2 * f3 * pix * piy - (f4 * pix * piy) / 2. + 2 * f6 * Power(pfz, 2) * pix * piy - 2 * f6 * pfy * pfz * pix * piz - 2 * f6 * pfx * pfz * piy * piz +
                  2 * f6 * pfx * pfy * Power(piz, 2);
            return cdouble(rep, imp);
        case 4:
            rep = -(f3 * pfx * pfz) - (f4 * pfx * pfz) / 4. + f3 * pfz * pix - (f4 * pfz * pix) / 4. + (f5 * pfz * pix) / 2. - f6 * pfy * pfz * pix * piy + f6 * pfx * pfz * Power(piy, 2) + f3 * pfx * piz - (f4 * pfx * piz) / 4. - (f5 * pfx * piz) / 2. - f3 * pix * piz - (f4 * pix * piz) / 4. +
                  f6 * Power(pfy, 2) * pix * piz - f6 * pfx * pfy * piy * piz;
            imp = -(f3 * pfy * pfz) - (f4 * pfy * pfz) / 4. + f6 * pfy * pfz * Power(pix, 2) + f3 * pfz * piy - (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. - f6 * pfx * pfz * pix * piy + f3 * pfy * piz - (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. - f6 * pfx * pfy * pix * piz - f3 * piy * piz -
                  (f4 * piy * piz) / 4. + f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 5:
            rep = f1 - f2 - f3 * Power(pfz, 2) - (f4 * Power(pfz, 2)) / 4. - f6 * Power(pfy, 2) * Power(pix, 2) + 2 * f6 * pfx * pfy * pix * piy - f6 * Power(pfx, 2) * Power(piy, 2) + 2 * f3 * pfz * piz - (f4 * pfz * piz) / 2. - f3 * Power(piz, 2) - (f4 * Power(piz, 2)) / 4.;
            imp = 0;
            return cdouble(rep, imp);
        case 6:
            rep = 2 * f2 + f3 * Power(pfx, 2) + (f4 * Power(pfx, 2)) / 4. + f3 * Power(pfy, 2) + (f4 * Power(pfy, 2)) / 4. - 2 * f3 * pfx * pix + (f4 * pfx * pix) / 2. + f3 * Power(pix, 2) + (f4 * Power(pix, 2)) / 4. + f6 * Power(pfz, 2) * Power(pix, 2) - 2 * f3 * pfy * piy + (f4 * pfy * piy) / 2. +
                  f3 * Power(piy, 2) + (f4 * Power(piy, 2)) / 4. + f6 * Power(pfz, 2) * Power(piy, 2) - 2 * f6 * pfx * pfz * pix * piz - 2 * f6 * pfy * pfz * piy * piz + f6 * Power(pfx, 2) * Power(piz, 2) + f6 * Power(pfy, 2) * Power(piz, 2);
            imp = 0;
            return cdouble(rep, imp);
        case 7:
            rep = f3 * pfx * pfz + (f4 * pfx * pfz) / 4. - f3 * pfz * pix + (f4 * pfz * pix) / 4. - (f5 * pfz * pix) / 2. + f6 * pfy * pfz * pix * piy - f6 * pfx * pfz * Power(piy, 2) - f3 * pfx * piz + (f4 * pfx * piz) / 4. + (f5 * pfx * piz) / 2. + f3 * pix * piz + (f4 * pix * piz) / 4. -
                  f6 * Power(pfy, 2) * pix * piz + f6 * pfx * pfy * piy * piz;
            imp = -(f3 * pfy * pfz) - (f4 * pfy * pfz) / 4. + f6 * pfy * pfz * Power(pix, 2) + f3 * pfz * piy - (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. - f6 * pfx * pfz * pix * piy + f3 * pfy * piz - (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. - f6 * pfx * pfy * pix * piz - f3 * piy * piz -
                  (f4 * piy * piz) / 4. + f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 8:
            rep = -(f3 * pfx * pfz) - (f4 * pfx * pfz) / 4. + f3 * pfz * pix - (f4 * pfz * pix) / 4. + (f5 * pfz * pix) / 2. - f6 * pfy * pfz * pix * piy + f6 * pfx * pfz * Power(piy, 2) + f3 * pfx * piz - (f4 * pfx * piz) / 4. - (f5 * pfx * piz) / 2. - f3 * pix * piz - (f4 * pix * piz) / 4. +
                  f6 * Power(pfy, 2) * pix * piz - f6 * pfx * pfy * piy * piz;
            imp = -(f3 * pfy * pfz) - (f4 * pfy * pfz) / 4. + f6 * pfy * pfz * Power(pix, 2) + f3 * pfz * piy - (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. - f6 * pfx * pfz * pix * piy + f3 * pfy * piz - (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. - f6 * pfx * pfy * pix * piz - f3 * piy * piz -
                  (f4 * piy * piz) / 4. + f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 9:
            rep = 2 * f2 + f3 * Power(pfx, 2) + (f4 * Power(pfx, 2)) / 4. + f3 * Power(pfy, 2) + (f4 * Power(pfy, 2)) / 4. - 2 * f3 * pfx * pix + (f4 * pfx * pix) / 2. + f3 * Power(pix, 2) + (f4 * Power(pix, 2)) / 4. + f6 * Power(pfz, 2) * Power(pix, 2) - 2 * f3 * pfy * piy + (f4 * pfy * piy) / 2. +
                  f3 * Power(piy, 2) + (f4 * Power(piy, 2)) / 4. + f6 * Power(pfz, 2) * Power(piy, 2) - 2 * f6 * pfx * pfz * pix * piz - 2 * f6 * pfy * pfz * piy * piz + f6 * Power(pfx, 2) * Power(piz, 2) + f6 * Power(pfy, 2) * Power(piz, 2);
            imp = 0;
            return cdouble(rep, imp);
        case 10:
            rep = f1 - f2 - f3 * Power(pfz, 2) - (f4 * Power(pfz, 2)) / 4. - f6 * Power(pfy, 2) * Power(pix, 2) + 2 * f6 * pfx * pfy * pix * piy - f6 * Power(pfx, 2) * Power(piy, 2) + 2 * f3 * pfz * piz - (f4 * pfz * piz) / 2. - f3 * Power(piz, 2) - (f4 * Power(piz, 2)) / 4.;
            imp = 0;
            return cdouble(rep, imp);
        case 11:
            rep = f3 * pfx * pfz + (f4 * pfx * pfz) / 4. - f3 * pfz * pix + (f4 * pfz * pix) / 4. - (f5 * pfz * pix) / 2. + f6 * pfy * pfz * pix * piy - f6 * pfx * pfz * Power(piy, 2) - f3 * pfx * piz + (f4 * pfx * piz) / 4. + (f5 * pfx * piz) / 2. + f3 * pix * piz + (f4 * pix * piz) / 4. -
                  f6 * Power(pfy, 2) * pix * piz + f6 * pfx * pfy * piy * piz;
            imp = -(f3 * pfy * pfz) - (f4 * pfy * pfz) / 4. + f6 * pfy * pfz * Power(pix, 2) + f3 * pfz * piy - (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. - f6 * pfx * pfz * pix * piy + f3 * pfy * piz - (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. - f6 * pfx * pfy * pix * piz - f3 * piy * piz -
                  (f4 * piy * piz) / 4. + f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 12:
            rep = f3 * Power(pfx, 2) + (f4 * Power(pfx, 2)) / 4. - f3 * Power(pfy, 2) - (f4 * Power(pfy, 2)) / 4. - 2 * f3 * pfx * pix + (f4 * pfx * pix) / 2. + f3 * Power(pix, 2) + (f4 * Power(pix, 2)) / 4. - f6 * Power(pfz, 2) * Power(pix, 2) + 2 * f3 * pfy * piy - (f4 * pfy * piy) / 2. - f3 * Power(piy, 2) -
                  (f4 * Power(piy, 2)) / 4. + f6 * Power(pfz, 2) * Power(piy, 2) + 2 * f6 * pfx * pfz * pix * piz - 2 * f6 * pfy * pfz * piy * piz - f6 * Power(pfx, 2) * Power(piz, 2) + f6 * Power(pfy, 2) * Power(piz, 2);
            imp = 2 * f3 * pfx * pfy + (f4 * pfx * pfy) / 2. - 2 * f3 * pfy * pix + (f4 * pfy * pix) / 2. - 2 * f3 * pfx * piy + (f4 * pfx * piy) / 2. + 2 * f3 * pix * piy + (f4 * pix * piy) / 2. - 2 * f6 * Power(pfz, 2) * pix * piy + 2 * f6 * pfy * pfz * pix * piz + 2 * f6 * pfx * pfz * piy * piz -
                  2 * f6 * pfx * pfy * Power(piz, 2);
            return cdouble(rep, imp);
        case 13:
            rep = f3 * pfx * pfz + (f4 * pfx * pfz) / 4. - f3 * pfz * pix + (f4 * pfz * pix) / 4. + (f5 * pfz * pix) / 2. + f6 * pfy * pfz * pix * piy - f6 * pfx * pfz * Power(piy, 2) - f3 * pfx * piz + (f4 * pfx * piz) / 4. - (f5 * pfx * piz) / 2. + f3 * pix * piz + (f4 * pix * piz) / 4. -
                  f6 * Power(pfy, 2) * pix * piz + f6 * pfx * pfy * piy * piz;
            imp = f3 * pfy * pfz + (f4 * pfy * pfz) / 4. - f6 * pfy * pfz * Power(pix, 2) - f3 * pfz * piy + (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. + f6 * pfx * pfz * pix * piy - f3 * pfy * piz + (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. + f6 * pfx * pfy * pix * piz + f3 * piy * piz + (f4 * piy * piz) / 4. -
                  f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 14:
            rep = f3 * pfx * pfz + (f4 * pfx * pfz) / 4. - f3 * pfz * pix + (f4 * pfz * pix) / 4. + (f5 * pfz * pix) / 2. + f6 * pfy * pfz * pix * piy - f6 * pfx * pfz * Power(piy, 2) - f3 * pfx * piz + (f4 * pfx * piz) / 4. - (f5 * pfx * piz) / 2. + f3 * pix * piz + (f4 * pix * piz) / 4. -
                  f6 * Power(pfy, 2) * pix * piz + f6 * pfx * pfy * piy * piz;
            imp = f3 * pfy * pfz + (f4 * pfy * pfz) / 4. - f6 * pfy * pfz * Power(pix, 2) - f3 * pfz * piy + (f4 * pfz * piy) / 4. + (f5 * pfz * piy) / 2. + f6 * pfx * pfz * pix * piy - f3 * pfy * piz + (f4 * pfy * piz) / 4. - (f5 * pfy * piz) / 2. + f6 * pfx * pfy * pix * piz + f3 * piy * piz + (f4 * piy * piz) / 4. -
                  f6 * Power(pfx, 2) * piy * piz;
            return cdouble(rep, imp);
        case 15:
            rep = f1 + f2 + f3 * Power(pfz, 2) + (f4 * Power(pfz, 2)) / 4. + f6 * Power(pfy, 2) * Power(pix, 2) - 2 * f6 * pfx * pfy * pix * piy + f6 * Power(pfx, 2) * Power(piy, 2) - 2 * f3 * pfz * piz + (f4 * pfz * piz) / 2. + f3 * Power(piz, 2) + (f4 * Power(piz, 2)) / 4.;
            imp = f5 * pfy * pix - f5 * pfx * piy;
            return cdouble(rep, imp);
        default: std::cerr << "Error in potential_auto!!!" << std::endl; std::exit(-1);
    }
}

} // namespace apwp

#endif // APWP_HPP