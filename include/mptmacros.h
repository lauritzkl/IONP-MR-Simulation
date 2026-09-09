/**
 * @file mptmacros.h
 *
 * Compiled physical constants, arithmetic helpers and explicit-cast macros.
 * Lengths are in metres, times in seconds and fields in tesla. These values
 * are part of the implemented model; changing them changes simulation results.
 */

/*****************************************************************************
 *
 *
 * @file   mptmacros.h
 * @author
 * @date
 *
 *****************************************************************************/
#ifndef MPTMACROS_H_
#  define MPTMACROS_H_ 1
/* ------------------------------------------------------------------------- */
const double D_PI = 3.14159265;
const double TWO_PI = D_PI * 2.0;
const double FOURTHIRD_PI = D_PI * 4.0 / 3.0;
const double FOUR_PI = D_PI * 4.0;

const double D_VOL_FRAC = D_PI*1.0E-6; // Default particle volume fraction.

const double D_D_CONST = 3E-9; // Bulk-water diffusion coefficient(m^2/s).

const double D_GI = 2.6752218744E8;// Proton gyromagnetic ratio, rad/(s*T); legacy NIST value.

//calculation of constant part of magnetic field
const double D_MO = FOUR_PI * 1.0E-7; // Permeability constant, T*m/A; legacy model value.

const double D_M_IONP = 3.8E5; // IONP magnetization(A/m); retained model parameter.

const double B_fixed = D_MO*D_M_IONP*(1./3.); // Dipole prefactor(T); multiplied by core radius cubed.

const double D_T_STEP = 10E-6;// Legacy step constant: 10 microseconds; not the configured runtime default.

template <class T> //instead of overload
T square(T a) {
    return a * a;
}

template <class T>
T cube(T a) {
    return a * a * a;
}
/* ------------------------------------------------------------------------- */
/** @name CAST macros                                                        */
/**@{*/
/** Macro to static cast a value to double type                              */
#  define SC_D( VAL ) ( static_cast< double >( VAL ) )
/** Macro to static cast a value to a long type                              */
#  define SC_L( VAL ) ( static_cast< long >( VAL ) )
/** Macro to static cast a value to a int32_t type                           */
#  define SC_I32( VAL ) ( static_cast< int32_t >( VAL ) )
/** Macro to static cast a value to a int16_t type                           */
#  define SC_I16( VAL ) ( static_cast< int16_t >( VAL ) )
/** Macro to static cast a value to a unsigned short type                    */
#  define SC_US( VAL ) ( static_cast< unsigned short >( VAL ) )
/** Macro to static cast a value to a unsigned int type                      */
#  define SC_UI( VAL ) ( static_cast< unsigned int >( VAL ) )
/** Macro to static cast a value to a int type                               */
#  define SC_I( VAL ) ( static_cast< int >( VAL ) )
/** Macro to static cast a value to a unsigned long(int) type               */
#  define SC_UL( VAL ) ( static_cast< unsigned long >( VAL ) )
/** Macro to static cast a value to a std::string type                       */
#  define SC_S( VAL ) ( static_cast< std::string >( VAL ) )
/**@}*/

/* ------------------------------------------------------------------------- */
#endif // #ifndef MPTMACROS_H_
/* ---- EOF ---------------------------------------------------------------- */
