using System;
using System.Text;

namespace Solar.Business
{
    /*============================================================================
    *    Contains:
    *        S_solpos     (computes solar position and intensity
    *                      from time and place)
    *
    *            INPUTS:     (via posdata struct) year, daynum, hour,
    *                        minute, second, latitude, longitude, timezone,
    *                        intervl
    *            OPTIONAL:   (via posdata struct) month, day, press, temp, tilt,
    *                        aspect, function
    *            OUTPUTS:    EVERY variable in the struct posdata
    *                            (defined in solpos.h)
    *
    *                       NOTE: Certain conditions exist during which some of
    *                       the output variables are undefined or cannot be
    *                       calculated.  In these cases, the variables are
    *                       returned with flag values indicating such.  In other
    *                       cases, the variables may return a realistic, though
    *                       invalid, value. These variables and the flag values
    *                       or invalid conditions are listed below:
    *
    *                       amass     -1.0 at zenetr angles greater than 93.0
    *                                 degrees
    *                       ampress   -1.0 at zenetr angles greater than 93.0
    *                                 degrees
    *                       azim      invalid at zenetr angle 0.0 or latitude
    *                                 +/-90.0 or at night
    *                       elevetr   limited to -9 degrees at night
    *                       etr       0.0 at night
    *                       etrn      0.0 at night
    *                       etrtilt   0.0 when cosinc is less than 0
    *                       prime     invalid at zenetr angles greater than 93.0
    *                                 degrees
    *                       sretr     +/- 2999.0 during periods of 24 hour sunup or
    *                                 sundown
    *                       ssetr     +/- 2999.0 during periods of 24 hour sunup or
    *                                 sundown
    *                       ssha      invalid at the North and South Poles
    *                       unprime   invalid at zenetr angles greater than 93.0
    *                                 degrees
    *                       zenetr    limited to 99.0 degrees at night
    *
    *        S_init       (optional initialization for all input parameters in
    *                      the posdata struct)
    *           INPUTS:     struct posdata*
    *           OUTPUTS:    struct posdata*
    *
    *                     (Note: initializes the required S_solpos INPUTS above
    *                      to out-of-bounds conditions, forcing the user to
    *                      supply the parameters; initializes the OPTIONAL
    *                      S_solpos inputs above to nominal values.)
    *
    *       S_decode      (optional utility for decoding the S_solpos return code)
    *           INPUTS:     long integer S_solpos return value, struct posdata*
    *           OUTPUTS:    text to stderr
    *
    *    Usage:
    *         In calling program, just after other 'includes', insert:
    *
    *              #include "solpos00.h"
    *
    *         Function calls:
    *              S_init(struct posdata*)  [optional]
    *              .
    *              .
    *              [set time and location parameters before S_solpos call]
    *              .
    *              .
    *              int retval = S_solpos(struct posdata*)
    *              S_decode(int retval, struct posdata*) [optional]
    *                  (Note: you should always look at the S_solpos return
    *                   value, which contains error codes. S_decode is one option
    *                   for examining these codes.  It can also serve as a
    *                   template for building your own application-specific
    *                   decoder.)
    *
    *    Martin Rymes
    *    National Renewable Energy Laboratory
    *    25 March 1998
    *
    *    27 April 1999 REVISION:  Corrected leap year in S_date.
    *    13 January 2000 REVISION:  SMW converted to structure posdata parameter
    *                               and subdivided into functions.
    *    01 February 2001 REVISION: SMW corrected ecobli calculation 
    *                               (changed sign). Error is small (max 0.015 deg
    *                               in calculation of declination angle)
    *----------------------------------------------------------------------------*/

    public class Solpos
    {

        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        *
        * Temporary global variables used only in this file:
        *
        *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        static int[,] month_days = { { 0,   0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334 },
                                     { 0,   0,  31,  60,  91, 121, 152, 182, 213, 244, 274, 305, 335 } };
        /* cumulative number of days prior to beginning of month */

        static double degrad = 57.295779513; /* converts from radians to degrees */
        static double raddeg = 0.0174532925; /* converts from degrees to radians */

        /*============================================================================
        *    Long integer function ExSolpos, adapted from the VAX solar libraries
        *
        *    This function calculates the apparent solar position and the
        *    intensity of the sun (theoretical maximum solar energy) from
        *    time and place on Earth.
        *
        *    Requires (from the struct posdata parameter):
        *        Date and time:
        *            year
        *            daynum   (requirement depends on the S_DOY switch)
        *            month    (requirement depends on the S_DOY switch)
        *            day      (requirement depends on the S_DOY switch)
        *            hour
        *            minute
        *            second
        *            interval  DEFAULT 0
        *        Location:
        *            latitude
        *            longitude
        *        Location/time adjuster:
        *            timezone
        *        Atmospheric pressure and temperature:
        *            press     DEFAULT 1013.0 mb
        *            temp      DEFAULT 10.0 degrees C
        *        Tilt of flat surface that receives solar energy:
        *            aspect    DEFAULT 180 (South)
        *            tilt      DEFAULT 0 (Horizontal)
        *        Function Switch (codes defined in solpos.h)
        *            function  DEFAULT S_ALL
        *
        *    Returns (via the struct posdata parameter):
        *        everything defined in the struct posdata in solpos.h.
        *----------------------------------------------------------------------------*/
        public ErrorCode ExSolpos(ref Postdata pdat)
        {
            ErrorCode retval;

            Trigdata tdat = new Trigdata();

            if ((retval = Validate(ref pdat)) != 0) // validate the inputs
                return retval;

            if (pdat.Function.HasFlag(Functions.L_DOY))
                Doy2Dom(ref pdat); // convert input doy to month-day
            else
                Dom2Doy(ref pdat); // convert input month-day to doy

            if (pdat.Function.HasFlag(Functions.L_GEOM))
                Geometry(ref pdat); // do basic geometry calculations

            if (pdat.Function.HasFlag(Functions.L_ZENETR)) // etr at non-refracted zenith angle
                Zen_no_ref(ref pdat, ref tdat);

            if (pdat.Function.HasFlag(Functions.L_SSHA)) // Sunset hour calculation
                Ssha(ref pdat, ref tdat);

            if (pdat.Function.HasFlag(Functions.L_SBCF)) // Shadowband correction factor
                Sbcf(ref pdat, ref tdat);

            if (pdat.Function.HasFlag(Functions.L_TST)) // true solar time
                Tst(ref pdat);

            if (pdat.Function.HasFlag(Functions.L_SRSS)) // sunrise/sunset calculations
                Srss(ref pdat);

            if (pdat.Function.HasFlag(Functions.L_SOLAZM)) // solar azimuth calculations
                Sazm(ref pdat, ref tdat);

            if (pdat.Function.HasFlag(Functions.L_REFRAC)) // atmospheric refraction calculations
                Refrac(ref pdat);

            if (pdat.Function.HasFlag(Functions.L_AMASS)) // airmass calculations
                Amass(ref pdat);

            if (pdat.Function.HasFlag(Functions.L_PRIME)) // kt-prime/unprime calculations
                Prime(ref pdat);

            if (pdat.Function.HasFlag(Functions.L_ETR)) // ETR and ETRN (refracted)
                Etr(ref pdat);

            if (pdat.Function.HasFlag(Functions.L_TILT)) // tilt calculations
                Tilt(ref pdat);

            return 0;
        }



        /*============================================================================
        *    Local long int function Validate
        *
        *    Validates the input parameters
        *----------------------------------------------------------------------------*/
        private static ErrorCode Validate(ref Postdata pdat)
        {
            ErrorCode retval = 0; // start with no errors

            /* No absurd dates, please. */
            if (pdat.Function.HasFlag(Functions.L_GEOM))
            {
                if (pdat.Year < 1950 || pdat.Year > 2050) // limits of algorithm
                    retval |= ErrorCode.S_YEAR_ERROR;

                if (!pdat.Function.HasFlag(Functions.S_DOY) && (pdat.Month < 1 || pdat.Month > 12))
                    retval |= ErrorCode.S_MONTH_ERROR;

                if (!pdat.Function.HasFlag(Functions.S_DOY) && (pdat.Day < 1 || pdat.Day > 31))
                    retval |= ErrorCode.S_DAY_ERROR;

                if (pdat.Function.HasFlag(Functions.S_DOY) && (pdat.Daynum < 1 || pdat.Daynum > 366))
                    retval |= ErrorCode.S_DOY_ERROR;


                /* No absurd times, please. */
                if (pdat.Hour < 0 || pdat.Hour > 24)
                    retval |= ErrorCode.S_HOUR_ERROR;

                if (pdat.Minute < 0 || pdat.Minute > 59)
                    retval |= ErrorCode.S_MINUTE_ERROR;

                if (pdat.Second < 0 || pdat.Second > 59)
                    retval |= ErrorCode.S_SECOND_ERROR;

                if (pdat.Hour == 24 && pdat.Minute > 0) // no more than 24 hrs
                    retval |= ErrorCode.S_HOUR_ERROR | ErrorCode.S_MINUTE_ERROR;

                if (pdat.Hour == 24 && pdat.Second > 0) // no more than 24 hrs
                    retval |= ErrorCode.S_HOUR_ERROR | ErrorCode.S_SECOND_ERROR;

                if (Math.Abs(pdat.Timezone) > 12.0)
                    retval |= ErrorCode.S_TZONE_ERROR;

                if (pdat.Interval < 0 || pdat.Interval > 28800)
                    retval |= ErrorCode.S_INTRVL_ERROR;


                /* No absurd locations, please. */
                if (Math.Abs(pdat.Longitude) > 180.0)
                    retval |= ErrorCode.S_LON_ERROR;

                if (Math.Abs(pdat.Latitude) > 90.0)
                    retval |= ErrorCode.S_LAT_ERROR;

            }

            /* No silly temperatures or pressures, please. */
            if (pdat.Function.HasFlag(Functions.L_REFRAC) && Math.Abs(pdat.Temp) > 100.0)
                retval |= ErrorCode.S_TEMP_ERROR;

            if (pdat.Function.HasFlag(Functions.L_REFRAC) && pdat.Press < 0.0 || pdat.Press > 2000.0)
                retval |= ErrorCode.S_PRESS_ERROR;


            /* No out of bounds tilts, please */
            if (pdat.Function.HasFlag(Functions.L_TILT) && Math.Abs(pdat.Tilt) > 180.0)
                retval |= ErrorCode.S_TILT_ERROR;

            if (pdat.Function.HasFlag(Functions.L_TILT) && Math.Abs(pdat.Aspect) > 360.0)
                retval |= ErrorCode.S_ASPECT_ERROR;


            /* No oddball shadowbands, please */
            if (pdat.Function.HasFlag(Functions.L_SBCF) && pdat.Sbwid < 1.0 || pdat.Sbwid > 100.0)
                retval |= ErrorCode.S_SBWID_ERROR;

            if (pdat.Function.HasFlag(Functions.L_SBCF) && pdat.Sbrad < 1.0 || pdat.Sbrad > 100.0)
                retval |= ErrorCode.S_SBRAD_ERROR;

            if (pdat.Function.HasFlag(Functions.L_SBCF) && Math.Abs(pdat.Sbsky) > 1.0)
                retval |= ErrorCode.S_SBSKY_ERROR;


            return retval;
        }


        /*============================================================================
        *    Local void function Doy2Dom
        *
        *    This function computes the month/day from the day number.
        *
        *    Requires (from struct posdata parameter):
        *        Year and day number:
        *            year
        *            daynum
        *
        *    Returns (via the struct posdata parameter):
        *            year
        *            month
        *            day
        *----------------------------------------------------------------------------*/
        private static void Doy2Dom(ref Postdata pdat)
        {
            int imon; // Month (month_days) array counter
            int leap; // leap year switch

            /* Set the leap year switch */
            if (((pdat.Year % 4) == 0) && (((pdat.Year % 100) != 0) || ((pdat.Year % 400) == 0)))
                leap = 1;
            else
                leap = 0;
            

            /* Find the month */
            imon = 12;
            while (pdat.Daynum <= month_days[leap, imon])
            {
                --imon;
            }

            /* Set the month and day of month */
            pdat.Month = imon;
            pdat.Day = pdat.Daynum - month_days[leap, imon];
        }

        /*============================================================================
        *    Local Void function Dom2Doy
        *
        *    Converts day-of-month to day-of-year
        *
        *    Requires (from struct posdata parameter):
        *            year
        *            month
        *            day
        *
        *    Returns (via the struct posdata parameter):
        *            year
        *            daynum
        *----------------------------------------------------------------------------*/
        private static void Dom2Doy(ref Postdata pdat)
        {
            pdat.Daynum = pdat.Day + month_days[0, pdat.Month];

            /* (adjust for leap year) */
            if (pdat.Year % 4 == 0 && (pdat.Year % 100 != 0 || pdat.Year % 400 == 0) && pdat.Month > 2)
                pdat.Daynum += 1;
            
        }

        /*============================================================================
        *    Local Void function Geometry
        *
        *    Does the underlying geometry for a given time and location
        *----------------------------------------------------------------------------*/
        private static void Geometry(ref Postdata pdat)
        {
            double bottom; // denominator (bottom) of the fraction
            double c2; // cosine of d2
            double cd; // cosine of the day angle or delination
            double d2; // pdat->dayang times two
            double delta; // difference between current year and 1949
            double s2; // sine of d2
            double sd; // sine of the day angle
            double top; // numerator (top) of the fraction
            int leap; // leap year counter

            /* Day angle */
            /*  Iqbal, M.  1983.  An Introduction to Solar Radiation.
                  Academic Press, NY., page 3 */
            pdat.Dayang = 360.0 * (pdat.Daynum - 1) / 365.0;

            /* Earth radius vector * solar constant = solar energy */
            /*  Spencer, J. W.  1971.  Fourier series representation of the
		        position of the sun.  Search 2 (5), page 172 */
            sd = Math.Sin(raddeg * pdat.Dayang);
            cd = Math.Cos(raddeg * pdat.Dayang);
            d2 = 2.0 * pdat.Dayang;
            c2 = Math.Cos(raddeg * d2);
            s2 = Math.Sin(raddeg * d2);

            pdat.Erv = 1.000110 + 0.034221 * cd + 0.001280 * sd;
            pdat.Erv += 0.000719 * c2 + 0.000077 * s2;

            /* Universal Coordinated (Greenwich standard) time */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Utime = pdat.Hour * 3600.0 + pdat.Minute * 60.0 + pdat.Second - pdat.Interval / 2.0;
            pdat.Utime = pdat.Utime / 3600.0 - pdat.Timezone;

            /* Julian Day minus 2,400,000 days (to eliminate roundoff errors) */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */

            /* No adjustment for century non-leap years since this function is
               bounded by 1950 - 2050 */
            delta = pdat.Year - 1949;
            leap = (int)(delta / 4.0);
            pdat.Julday = 32916.5 + delta * 365.0 + leap + pdat.Daynum + pdat.Utime / 24.0;

            /* Time used in the calculation of ecliptic coordinates */
            /* Noon 1 JAN 2000 = 2,400,000 + 51,545 days Julian Date */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Ectime = pdat.Julday - 51545.0;

            /* Mean longitude */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Mnlong = 280.460 + 0.9856474 * pdat.Ectime;

            /* (dump the multiples of 360, so the answer is between 0 and 360) */
            pdat.Mnlong -= 360.0 * (int)(pdat.Mnlong / 360.0);
            if (pdat.Mnlong < 0.0)
                pdat.Mnlong += 360.0;           

            /* Mean anomaly */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Mnanom = 357.528 + 0.9856003 * pdat.Ectime;

            /* (dump the multiples of 360, so the answer is between 0 and 360) */
            pdat.Mnanom -= 360.0 * (int)(pdat.Mnanom / 360.0);
            if (pdat.Mnanom < 0.0)
                pdat.Mnanom += 360.0;           

            /* Ecliptic longitude */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Eclong = pdat.Mnlong + 1.915 * Math.Sin(pdat.Mnanom * raddeg) + 0.020 * Math.Sin(2.0 * pdat.Mnanom * raddeg);

            /* (dump the multiples of 360, so the answer is between 0 and 360) */
            pdat.Eclong -= 360.0 * (int)(pdat.Eclong / 360.0);
            if (pdat.Eclong < 0.0)
                pdat.Eclong += 360.0;            

            /* Obliquity of the ecliptic */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */

            /* 02 Feb 2001 SMW corrected sign in the following line */
            /*  pdat->ecobli = 23.439 + 4.0e-07 * pdat->ectime;     */
            pdat.Ecobli = 23.439 - 4.0e-07 * pdat.Ectime;

            /* Declination */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Declin = degrad * Math.Asin(Math.Sin(pdat.Ecobli * raddeg) * Math.Sin(pdat.Eclong * raddeg));

            /* Right ascension */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            top = Math.Cos(raddeg * pdat.Ecobli) * Math.Sin(raddeg * pdat.Eclong);
            bottom = Math.Cos(raddeg * pdat.Eclong);

            pdat.Rascen = degrad * Math.Atan2(top, bottom);

            /* (make it a positive angle) */
            if (pdat.Rascen < 0.0)
                pdat.Rascen += 360.0;


            /* Greenwich mean sidereal time */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Gmst = 6.697375 + 0.0657098242 * pdat.Ectime + pdat.Utime;

            /* (dump the multiples of 24, so the answer is between 0 and 24) */
            pdat.Gmst -= 24.0 * (int)(pdat.Gmst / 24.0);
            if (pdat.Gmst < 0.0)
                pdat.Gmst += 24.0;


            /* Local mean sidereal time */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Lmst = pdat.Gmst * 15.0 + pdat.Longitude;

            /* (dump the multiples of 360, so the answer is between 0 and 360) */
            pdat.Lmst -= 360.0 * (int)(pdat.Lmst / 360.0);
            if (pdat.Lmst < 0.0)
                pdat.Lmst += 360.0;
            

            /* Hour angle */
            /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
		        approximate solar position (1950-2050).  Solar Energy 40 (3),
		        pp. 227-235. */
            pdat.Hrang = pdat.Lmst - pdat.Rascen;

            /* (force it between -180 and 180 degrees) */
            if (pdat.Hrang < -180.0)
                pdat.Hrang += 360.0;
            else if (pdat.Hrang > 180.0)
                pdat.Hrang -= 360.0;            
        }

        /*============================================================================
        *    Local Void function Zen_no_ref
        *
        *    ETR solar zenith angle
        *       Iqbal, M.  1983.  An Introduction to Solar Radiation.
        *            Academic Press, NY., page 15
        *----------------------------------------------------------------------------*/
        private static void Zen_no_ref(ref Postdata pdat, ref Trigdata tdat)
        {
            double cz; // cosine of the solar zenith angle

            Localtrig(ref pdat, ref tdat);
            cz = tdat.Sd * tdat.Sl + tdat.Cd * tdat.Cl * tdat.Ch;

            /* (watch out for the roundoff errors) */
            if (Math.Abs(cz) > 1.0)
            {
                if (cz >= 0.0)
                    cz = 1.0;
                else
                    cz = -1.0;
            }

            pdat.Zenetr = Math.Acos(cz) * degrad;

            /* (limit the degrees below the horizon to 9 [+90 -> 99]) */
            if (pdat.Zenetr > 99.0)
                pdat.Zenetr = 99.0;
            

            pdat.Elevetr = 90.0 - pdat.Zenetr;
        }


        /*============================================================================
        *    Local Void function Localtrig
        *
        *    Does trig on internal variable used by several functions
        *----------------------------------------------------------------------------*/
        private static void Localtrig(ref Postdata pdat, ref Trigdata tdat)
        {

            if (tdat.Sd < -900.0) // sd was initialized -999 as flag
            {
                tdat.Sd = 1.0; // reflag as having completed calculations
                if ((pdat.Function | Functions.CD_MASK) != 0)
                    tdat.Cd = Math.Cos(raddeg * pdat.Declin);

                if ((pdat.Function | Functions.CH_MASK) != 0)
                    tdat.Ch = Math.Cos(raddeg * pdat.Hrang);

                if ((pdat.Function | Functions.CL_MASK) != 0)
                    tdat.Cl = Math.Cos(raddeg * pdat.Latitude);

                if ((pdat.Function | Functions.SD_MASK) != 0)
                    tdat.Sd = Math.Sin(raddeg * pdat.Declin);

                if ((pdat.Function | Functions.SL_MASK) != 0)
                    tdat.Sl = Math.Sin(raddeg * pdat.Latitude);                
            }
        }

        /*============================================================================
        *    Local Void function Ssha
        *
        *    Sunset hour angle, degrees
        *       Iqbal, M.  1983.  An Introduction to Solar Radiation.
        *            Academic Press, NY., page 16
        *----------------------------------------------------------------------------*/
        private static void Ssha(ref Postdata pdat, ref Trigdata tdat)
        {
            double cssha; // cosine of the sunset hour angle
            double cdcl; // ( cd * cl )

            Localtrig(ref pdat, ref tdat);
            cdcl = tdat.Cd * tdat.Cl;

            if (Math.Abs(cdcl) >= 0.001)
            {
                cssha = -tdat.Sl * tdat.Sd / cdcl;

                /* This keeps the cosine from blowing on roundoff */
                if (cssha < -1.0)
                    pdat.Ssha = 180.0;
                else if (cssha > 1.0)
                    pdat.Ssha = 0.0;
                else
                    pdat.Ssha = degrad * Math.Acos(cssha);
            }
            else if (((pdat.Declin >= 0.0) && (pdat.Latitude > 0.0)) || ((pdat.Declin < 0.0) && (pdat.Latitude < 0.0)))
                pdat.Ssha = 180.0;
            else
                pdat.Ssha = 0.0;
        }

        /*============================================================================
        *    Local Void function Sbcf
        *
        *    Shadowband correction factor
        *       Drummond, A. J.  1956.  A contribution to absolute pyrheliometry.
        *            Q. J. R. Meteorol. Soc. 82, pp. 481-493
        *----------------------------------------------------------------------------*/
        private static void Sbcf(ref Postdata pdat, ref Trigdata tdat)
        {
            double p; // used to compute sbcf
            double t1;
            double t2;

            Localtrig(ref pdat, ref tdat);
            p = 0.6366198 * pdat.Sbwid / pdat.Sbrad * Math.Pow(tdat.Cd, 3);
            t1 = tdat.Sl * tdat.Sd * pdat.Ssha * raddeg;
            t2 = tdat.Cl * tdat.Cd * Math.Sin(pdat.Ssha * raddeg);
            pdat.Sbcf = pdat.Sbsky + 1.0 / (1.0 - p * (t1 + t2));
        }

        /*============================================================================
        *    Local Void function Tst
        *
        *    TST -> True Solar Time = local standard time + TSTfix, time
        *      in minutes from midnight.
        *        Iqbal, M.  1983.  An Introduction to Solar Radiation.
        *            Academic Press, NY., page 13
        *----------------------------------------------------------------------------*/
        private static void Tst(ref Postdata pdat)
        {
            pdat.Tst = (180.0 + pdat.Hrang) * 4.0;
            pdat.Tstfix = pdat.Tst - pdat.Hour * 60.0 - pdat.Minute - pdat.Second / 60.0 + pdat.Interval / 120.0; // add back half of the interval

            /* bound tstfix to this day */
            while (pdat.Tstfix > 720.0)
                pdat.Tstfix -= 1440.0;
            
            while (pdat.Tstfix < -720.0)
                pdat.Tstfix += 1440.0;
            
            pdat.Eqntim = pdat.Tstfix + 60.0 * pdat.Timezone - 4.0 * pdat.Longitude;
        }

        /*============================================================================
        *    Local Void function Srss
        *
        *    Sunrise and sunset times (minutes from midnight)
        *----------------------------------------------------------------------------*/
        private static void Srss(ref Postdata pdat)
        {
            if (pdat.Ssha <= 1.0)
            {
                pdat.Sretr = 2999.0;
                pdat.Ssetr = -2999.0;
            }
            else if (pdat.Ssha >= 179.0)
            {
                pdat.Sretr = -2999.0;
                pdat.Ssetr = 2999.0;
            }
            else
            {
                pdat.Sretr = 720.0 - 4.0 * pdat.Ssha - pdat.Tstfix;
                pdat.Ssetr = 720.0 + 4.0 * pdat.Ssha - pdat.Tstfix;
            }
        }

        /*============================================================================
        *    Local Void function Sazm
        *
        *    Solar azimuth angle
        *       Iqbal, M.  1983.  An Introduction to Solar Radiation.
        *            Academic Press, NY., page 15
        *----------------------------------------------------------------------------*/
        private static void Sazm(ref Postdata pdat, ref Trigdata tdat)
        {
            double ca; // cosine of the solar azimuth angle
            double ce; // cosine of the solar elevation
            double cecl; // ( ce * cl )
            double se; // sine of the solar elevation

            Localtrig(ref pdat, ref tdat);
            ce = Math.Cos(raddeg * pdat.Elevetr);
            se = Math.Sin(raddeg * pdat.Elevetr);

            pdat.Azim = 180.0;
            cecl = ce * tdat.Cl;
            if (Math.Abs(cecl) >= 0.001)
            {
                ca = (se * tdat.Sl - tdat.Sd) / cecl;
                if (ca > 1.0)
                {
                    ca = 1.0;
                }
                else if (ca < -1.0)
                {
                    ca = -1.0;
                }

                pdat.Azim = 180.0 - Math.Acos(ca) * degrad;
                if (pdat.Hrang > 0)
                {
                    pdat.Azim = 360.0 - pdat.Azim;
                }
            }
        }


        /*============================================================================
        *    Local Int function Refrac
        *
        *    Refraction correction, degrees
        *        Zimmerman, John C.  1981.  Sun-pointing programs and their
        *            accuracy.
        *            SAND81-0761, Experimental Systems Operation Division 4721,
        *            Sandia National Laboratories, Albuquerque, NM.
        *----------------------------------------------------------------------------*/
        private static void Refrac(ref Postdata pdat)
        {
            double prestemp; // temporary pressure/temperature correction
            double refcor; // temporary refraction correction
            double tanelev; // tangent of the solar elevation angle

            /* If the sun is near zenith, the algorithm bombs; refraction near 0 */
            if (pdat.Elevetr > 85.0)
                refcor = 0.0;           
            /* Otherwise, we have refraction */
            else
            {
                tanelev = Math.Tan(raddeg * pdat.Elevetr);
                if (pdat.Elevetr >= 5.0)
                {
                    refcor = 58.1 / tanelev - 0.07 / (Math.Pow(tanelev, 3)) + 0.000086 / (Math.Pow(tanelev, 5));
                }
                else if (pdat.Elevetr >= -0.575)
                {
                    refcor = 1735.0 + pdat.Elevetr * (-518.2 + pdat.Elevetr * (103.4 + pdat.Elevetr * (-12.79 + pdat.Elevetr * 0.711)));
                }
                else
                {
                    refcor = -20.774 / tanelev;
                }

                prestemp = (pdat.Press * 283.0) / (1013.0 * (273.0 + pdat.Temp));
                refcor *= prestemp / 3600.0;
            }

            /* Refracted solar elevation angle */
            pdat.Elevref = pdat.Elevetr + refcor;

            /* (limit the degrees below the horizon to 9) */
            if (pdat.Elevref < -9.0)
            {
                pdat.Elevref = -9.0;
            }

            /* Refracted solar zenith angle */
            pdat.Zenref = 90.0 - pdat.Elevref;
            pdat.Coszen = Math.Cos(raddeg * pdat.Zenref);
        }

        /*============================================================================
        *    Local Void function  Amass
        *
        *    Airmass
        *       Kasten, F. and Young, A.  1989.  Revised optical air mass
        *            tables and approximation formula.  Applied Optics 28 (22),
        *            pp. 4735-4738
        *----------------------------------------------------------------------------*/
        private static void Amass(ref Postdata pdat)
        {
            if (pdat.Zenref > 93.0)
            {
                pdat.Amass = -1.0;
                pdat.Ampress = -1.0;
            }
            else
            {
                pdat.Amass = 1.0 / (Math.Cos(raddeg * pdat.Zenref) + 0.50572 * Math.Pow((96.07995 - pdat.Zenref), -1.6364));

                pdat.Ampress = pdat.Amass * pdat.Press / 1013.0;
            }
        }


        /*============================================================================
        *    Local Void function Prime
        *
        *    Prime and Unprime
        *    Prime  converts Kt to normalized Kt', etc.
        *       Unprime deconverts Kt' to Kt, etc.
        *            Perez, R., P. Ineichen, Seals, R., & Zelenka, A.  1990.  Making
        *            full use of the clearness index for parameterizing hourly
        *            insolation conditions. Solar Energy 45 (2), pp. 111-114
        *----------------------------------------------------------------------------*/
        private static void Prime(ref Postdata pdat)
        {
            pdat.Unprime = 1.031 * Math.Exp(-1.4 / (0.9 + 9.4 / pdat.Amass)) + 0.1;
            pdat.Prime = 1.0 / pdat.Unprime;
        }


        /*============================================================================
        *    Local Void function Etr
        *
        *    Extraterrestrial (top-of-atmosphere) solar irradiance
        *----------------------------------------------------------------------------*/
        private static void Etr(ref Postdata pdat)
        {
            if (pdat.Coszen > 0.0)
            {
                pdat.Etrn = pdat.Solcon * pdat.Erv;
                pdat.Etr = pdat.Etrn * pdat.Coszen;
            }
            else
            {
                pdat.Etrn = 0.0;
                pdat.Etr = 0.0;
            }
        }

        /*============================================================================
        *    Local Void function Tilt
        *
        *    ETR on a tilted surface
        *----------------------------------------------------------------------------*/
        private static void Tilt(ref Postdata pdat)
        {
            double ca; // cosine of the solar azimuth angle
            double cp; // cosine of the panel aspect
            double ct; // cosine of the panel tilt
            double sa; // sine of the solar azimuth angle
            double sp; // sine of the panel aspect
            double st; // sine of the panel tilt
            double sz; // sine of the refraction corrected solar zenith angle


            /* Cosine of the angle between the sun and a tipped flat surface,
               useful for calculating solar energy on tilted surfaces */
            ca = Math.Cos(raddeg * pdat.Azim);
            cp = Math.Cos(raddeg * pdat.Aspect);
            ct = Math.Cos(raddeg * pdat.Tilt);
            sa = Math.Sin(raddeg * pdat.Azim);
            sp = Math.Sin(raddeg * pdat.Aspect);
            st = Math.Sin(raddeg * pdat.Tilt);
            sz = Math.Sin(raddeg * pdat.Zenref);
            pdat.Cosinc = pdat.Coszen * ct + sz * st * (ca * cp + sa * sp);

            if (pdat.Cosinc > 0.0)
            {
                pdat.Etrtilt = pdat.Etrn * pdat.Cosinc;
            }
            else
            {
                pdat.Etrtilt = 0.0;
            }

        }


        /*============================================================================
        *    Void function S_decode
        *
        *    This function decodes the error codes from S_solpos return value
        *
        *    Requires the long integer return value from S_solpos
        *
        *    Returns descriptive text to stderr
        *----------------------------------------------------------------------------*/
        public string S_decode(ErrorCode code, ref Postdata pdat)
        {
            StringBuilder sb = new StringBuilder();

            if (code.HasFlag(ErrorCode.S_YEAR_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the year: {pdat.Year:D} [1950-2050]\n");

            if (code.HasFlag(ErrorCode.S_MONTH_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the month: {pdat.Month:D}\n");

            if (code.HasFlag(ErrorCode.S_DAY_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the day-of-month: {pdat.Day:D}\n");
        
            if (code.HasFlag(ErrorCode.S_DOY_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the day-of-year: {pdat.Daynum:D}\n");

            if (code.HasFlag(ErrorCode.S_HOUR_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the hour: {pdat.Hour:D}\n");
            
            if (code.HasFlag(ErrorCode.S_MINUTE_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the minute: {pdat.Minute:D}\n");

            if (code.HasFlag(ErrorCode.S_SECOND_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the second: {pdat.Second:D}\n");

            if (code.HasFlag(ErrorCode.S_TZONE_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the time zone: {pdat.Timezone:f}\n");

            if (code.HasFlag(ErrorCode.S_INTRVL_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the interval: {pdat.Interval:D}\n");

            if (code.HasFlag(ErrorCode.S_LAT_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the latitude: {pdat.Latitude:f}\n");

            if (code.HasFlag(ErrorCode.S_LON_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the longitude: {pdat.Longitude:f}\n");

            if (code.HasFlag(ErrorCode.S_TEMP_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the temperature: {pdat.Temp:f}\n");

            if (code.HasFlag(ErrorCode.S_PRESS_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the pressure: {pdat.Press:f}\n");

            if (code.HasFlag(ErrorCode.S_TILT_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the tilt: {pdat.Tilt:f}\n");

            if (code.HasFlag(ErrorCode.S_ASPECT_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the aspect: {pdat.Aspect:f}\n");

            if (code.HasFlag(ErrorCode.S_SBWID_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the shadowband width: {pdat.Sbwid:f}\n");

            if (code.HasFlag(ErrorCode.S_SBRAD_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the shadowband radius: {pdat.Sbrad:f}\n");

            if (code.HasFlag(ErrorCode.S_SBSKY_ERROR))
                sb.AppendLine($"S_decode ==> Please fix the shadowband sky factor: {pdat.Sbsky:f}\n");

            return sb.ToString();
        }



    }


    /*============================================================================
    *
    *     Enumerate the error codes
    *     (Bit positions are from least significant to most significant)
    *
    *----------------------------------------------------------------------------*/

    [Flags]
    public enum ErrorCode
    {   
        /*    Code          Bit                     Parameter            Range
        =============       =======               ===================   =============   */
        S_YEAR_ERROR =      1 << 0,         //  0   year                  1950 -  2050
        S_MONTH_ERROR =     1 << 1,         //  1   month                    1 -    12
        S_DAY_ERROR =       1 << 2,         //  2   day-of-month             1 -    31
        S_DOY_ERROR =       1 << 3,         //  3   day-of-year              1 -   366
        S_HOUR_ERROR =      1 << 4,         //  4   hour                     0 -    24
        S_MINUTE_ERROR =    1 << 5,         //  5   minute                   0 -    59
        S_SECOND_ERROR =    1 << 6,         //  6   second                   0 -    59
        S_TZONE_ERROR =     1 << 7,         //  7   time zone              -12 -    12
        S_INTRVL_ERROR =    1 << 8,         //  8   interval (seconds)       0 - 28800
        S_LAT_ERROR =       1 << 9,         //  9   latitude               -90 -    90
        S_LON_ERROR =       1 << 10,        // 10   longitude             -180 -   180
        S_TEMP_ERROR =      1 << 11,        // 11   temperature (deg. C)  -100 -   100
        S_PRESS_ERROR =     1 << 12,        // 12   pressure (millibars)     0 -  2000
        S_TILT_ERROR =      1 << 13,        // 13   tilt                   -90 -    90
        S_ASPECT_ERROR =    1 << 14,        // 14   aspect                -360 -   360
        S_SBWID_ERROR =     1 << 15,        // 15   shadow band width (cm)   1 -   100
        S_SBRAD_ERROR =     1 << 16,        // 16   shadow band radius (cm)  1 -   100
        S_SBSKY_ERROR =     1 << 17         // 17   shadow band sky factor  -1 -     1
    }


    [Flags]
    public enum Functions
    {
        L_DOY = 1 << 1,     // L_DOY = 0x0001;
        L_GEOM = 1 << 2,    //L_GEOM = 0x0002;
        L_ZENETR = 1 << 3,  //L_ZENETR = 0x0004;
        L_SSHA = 1 << 4,    //L_SSHA = 0x0008;
        L_SBCF = 1 << 5,    //L_SBCF = 0x0010;
        L_TST = 1 << 6,     //L_TST = 0x0020;
        L_SRSS = 1 << 7,    // L_SRSS = 0x0040;
        L_SOLAZM = 1 << 8,  //L_SOLAZM = 0x0080;
        L_REFRAC = 1 << 9,  //L_REFRAC = 0x0100; 
        L_AMASS = 1 << 10,  //L_AMASS = 0x0200;
        L_PRIME = 1 << 11,  //L_PRIME = 0x0400;
        L_TILT = 1 << 12,   //L_TILT = 0x0800;
        L_ETR = 1 << 13,    //L_ETR = 0x1000;

        L_DEFAULT = L_GEOM | L_ZENETR | L_SSHA | L_SBCF | L_TST | L_SRSS | L_SOLAZM | L_REFRAC | L_AMASS | L_PRIME | L_TILT | L_ETR, 
        L_ALL = L_DOY | L_GEOM | L_ZENETR | L_SSHA | L_SBCF | L_TST | L_SRSS | L_SOLAZM | L_REFRAC | L_AMASS | L_PRIME | L_TILT | L_ETR, // L_ALL = 0xFFFF;
        
        S_DOY = L_DOY,
        S_GEOM = L_GEOM | S_DOY,
        S_ZENETR = L_ZENETR | S_GEOM,
        S_SSHA = L_SSHA | S_GEOM,
        S_SBCF = L_SBCF | S_SSHA,
        S_TST = L_TST | S_GEOM,
        S_SRSS = L_SRSS | S_SSHA | S_TST,
        S_SOLAZM = L_SOLAZM | S_ZENETR,
        S_REFRAC = L_REFRAC | S_ZENETR,
        S_AMASS = L_AMASS | S_REFRAC,
        S_PRIME = L_PRIME | S_AMASS,
        S_TILT = L_TILT | S_SOLAZM | S_REFRAC,
        S_ETR = L_ETR | S_REFRAC,
        S_ALL = L_ALL,

        SD_MASK = L_ZENETR | L_SSHA | S_SBCF | S_SOLAZM,
        SL_MASK = L_ZENETR | L_SSHA | S_SBCF | S_SOLAZM,
        CL_MASK = L_ZENETR | L_SSHA | S_SBCF | S_SOLAZM,
        CD_MASK = L_ZENETR | L_SSHA | S_SBCF,
        CH_MASK = L_ZENETR
    }


    public class Trigdata // used to pass calculated values locally
    {
        public Trigdata()
        {
            /* initialize the trig structure */
            Sd = -999.0; // flag to force calculation of trig data
            Cd = 1.0;
            Ch = 1.0; // set the rest of these to something safe
            Cl = 1.0;
            Sl = 1.0;
        }

        public double Cd { get; set; } // cosine of the declination
        public double Ch { get; set; } // cosine of the hour angle
        public double Cl { get; set; } // cosine of the latitude
        public double Sd { get; set; } // sine of the declination
        public double Sl { get; set; } // sine of the latitude
    }

    public class Postdata
    {
        /***** ALPHABETICAL LIST OF COMMON VARIABLES *****/
        /* Each comment begins with a 1-column letter code:
            I:  INPUT variable
            O:  OUTPUT variable
            T:  TRANSITIONAL variable used in the algorithm, of interest only to the solar radiation
                modelers, and available to you because you may be one of them.

            The FUNCTION column indicates which sub-function within solpos must be switched on using the
            "function" parameter to calculate the desired output variable.  All function codes are
            defined in the solpos.h file.  The default S_ALL switch calculates all output variables.
            Multiple functions may be or'd to create a composite function switch.  For example,
            (S_TST | S_SBCF). Specifying only the functions for required output variables may allow solpos
            to execute more quickly.

            The S_DOY mask works as a toggle between the input date represented as a day number (daynum)
            or as month and day.  To set the switch (to use daynum input), the function is or'd; to
            clear the switch (to use month and day input), the function is inverted and and'd.

            For example:
                pdat->function |= S_DOY (sets daynum input)
                pdat->function &= ~S_DOY (sets month and day input)

            Whichever date form is used, S_solpos will calculate and return the variables(s) of the
            other form.  See the soltest.c program for other examples. */

        /* VARIABLE        I/O          Function   Description */
        /* -------------  -----------  ----------  ---------------------------------------*/

        public Postdata()
        {
            Day = -99;                  /* Day of month (May 27 = 27, etc.) */
            Daynum = -999;              /* Day number (day of year; Feb 1 = 32 ) */
            Hour = -99;                 /* Hour of day, 0 - 23 */
            Minute = -99;               /* Minute of hour, 0 - 59 */
            Month = -99;                /* Month number (Jan = 1, Feb = 2, etc.) */
            Second = -99;               /* Second of minute, 0 - 59 */
            Year = -99;                 /* 4-digit year */
            Interval = 0;               /* instantaneous measurement interval */
            Aspect = 180.0;             /* Azimuth of panel surface (direction it faces) N=0, E=90, S=180, W=270 */
            Latitude = -99.0;           /* Latitude, degrees north (south negative) */
            Longitude = -999.0;         /* Longitude, degrees east (west negative) */
            Press = 1013.0;             /* Surface pressure, millibars */
            Solcon = 1367.0;            /* Solar constant, 1367 W/sq m */
            Temp = 15.0;                /* Ambient dry-bulb temperature, degrees C */
            Tilt = 0.0;                 /* Degrees tilt from horizontal of panel */
            Timezone = -99.0;           /* Time zone, east (west negative). */
            Sbwid = 7.6;                /* Eppley shadow band width */
            Sbrad = 31.7;               /* Eppley shadow band radius */
            Sbsky = 0.04;               /* Drummond factor for partly cloudy skies */
            Function = Functions.S_ALL; /* compute all parameters */
        }

        public int Day { get; set; } /* I/O:       S_DOY Day of month (May 27 = 27, etc.) solpos will CALCULATE this by default,
                                                   or will optionally require it as input depending on the setting of the S_DOY
                                                   function switch. */
        public int Daynum { get; set; } /* I/O:    S_DOY Day number (day of year; Feb 1 = 32 ) solpos REQUIRES this by default, but
                                                   will optionally calculate it from month and day depending on the setting
                                                   of the S_DOY function switch. */
        public Functions Function { get; set; } /* I: Switch to choose functions for desired output. */
        public int Hour { get; set; } // I:        Hour of day, 0 - 23, DEFAULT = 12
        public int Interval { get; set; } /* I:    Interval of a measurement period in seconds.  Forces solpos to use the
                                                   time and date from the interval midpoint. The INPUT time (hour,
                                                   minute, and second) is assumed to be the END of the measurement
                                                   interval. */
        public int Minute { get; set; } // I:      Minute of hour, 0 - 59, DEFAULT = 0
        public int Month { get; set; } /* I/O:     S_DOY Month number (Jan = 1, Feb = 2, etc.) solpos will CALCULATE this by default,
                                                   or will optionally require it as input depending on the setting of the S_DOY
                                                   function switch. */
        public int Second { get; set; } // I:      Second of minute, 0 - 59, DEFAULT = 0
        public int Year { get; set; } /* I:        4-digit year (2-digit year is NOT allowed */

        /***** double's *****/
        public double Amass { get; set; } // O:      S_AMASS    Relative optical airmass
        public double Ampress { get; set; } // O:    S_AMASS    Pressure-corrected airmass
        public double Aspect { get; set; } /* I:     Azimuth    of panel surface (direction it faces) N=0, E=90, S=180, W=270, DEFAULT = 180 */
        public double Azim { get; set; } /* O:       S_SOLAZM   Solar azimuth angle:  N=0, E=90, S=180, W=270 */
        public double Cosinc { get; set; } /* O:     S_TILT     Cosine of solar incidence angle on panel */
        public double Coszen { get; set; } /* O:     S_REFRAC   Cosine of refraction corrected solar zenith angle */
        public double Dayang { get; set; } /* T:     S_GEOM     Day angle (daynum*360/year-length) degrees */
        public double Declin { get; set; } /* T:     S_GEOM     Declination--zenith angle of solar noon at equator, degrees NORTH */
        public double Eclong { get; set; } // T:     S_GEOM     Ecliptic longitude, degrees
        public double Ecobli { get; set; } // T:     S_GEOM     Obliquity of ecliptic
        public double Ectime { get; set; } // T:     S_GEOM     Time of ecliptic calculations
        public double Elevetr { get; set; } /* O:    S_ZENETR   Solar elevation, no atmospheric correction (= ETR) */
        public double Elevref { get; set; } /* O:    S_REFRAC   Solar elevation angle, deg. from horizon, refracted */
        public double Eqntim { get; set; } // T:     S_TST      Equation of time (TST - LMT), minutes
        public double Erv { get; set; } /* T:        S_GEOM     Earth radius vector (multiplied to solar constant) */
        public double Etr { get; set; } /* O:        S_ETR      Extraterrestrial (top-of-atmosphere) W/sq m global horizontal solar irradiance */
        public double Etrn { get; set; } /* O:       S_ETR      Extraterrestrial (top-of-atmosphere) W/sq m direct normal solar  irradiance */
        public double Etrtilt { get; set; } /* O:    S_TILT     Extraterrestrial (top-of-atmosphere) W/sq m global irradiance on a tilted surface */
        public double Gmst { get; set; } // T:       S_GEOM     Greenwich mean sidereal time, hours
        public double Hrang { get; set; } /* T:      S_GEOM     Hour angle--hour of sun from solar noon, degrees WEST */
        public double Julday { get; set; } /* T:     S_GEOM     Julian Day of 1 JAN 2000 minus 2,400,000 days (in order to regain single precision) */
        public double Latitude { get; set; } // I:              Latitude, degrees north (south negative)
        public double Longitude { get; set; } // I:             Longitude, degrees east (west negative)
        public double Lmst { get; set; } // T:       S_GEOM     Local mean sidereal time, degrees
        public double Mnanom { get; set; } // T:     S_GEOM     Mean anomaly, degrees
        public double Mnlong { get; set; } // T:     S_GEOM     Mean longitude, degrees
        public double Rascen { get; set; } // T:     S_GEOM     Right ascension, degrees
        public double Press { get; set; } /* I:                 Surface pressure, millibars (hPa), used for refraction correction and ampress */
        public double Prime { get; set; } // O:      S_PRIME    Factor that normalizes Kt, Kn, etc.
        public double Sbcf { get; set; } // O:       S_SBCF     Shadow-band correction factor
        public double Sbwid { get; set; } // I:                 Shadow-band width (cm)
        public double Sbrad { get; set; } // I:                 Shadow-band radius (cm)
        public double Sbsky { get; set; } // I:                 Shadow-band sky factor
        public double Solcon { get; set; } // I:                Solar constant (NREL uses 1367 W/sq m)
        public double Ssha { get; set; } // T:       S_SRHA     Sunset(/rise) hour angle, degrees
        public double Sretr { get; set; } /* O:      S_SRSS     Sunrise time, minutes from midnight, local, WITHOUT refraction */
        public double Ssetr { get; set; } /* O:      S_SRSS     Sunset time, minutes from midnight, local, WITHOUT refraction */
        public double Temp { get; set; } /* I:                  Ambient dry-bulb temperature, degrees C, used for refraction correction */
        public double Tilt { get; set; } // I:                  Degrees tilt from horizontal of panel
        public double Timezone { get; set; } /* I:              Time zone, east (west negative). USA:  Mountain = -7, Central = -6, etc. */
        public double Tst { get; set; } // T:        S_TST      True solar time, minutes from midnight
        public double Tstfix { get; set; } // T:     S_TST      True solar time - local standard time
        public double Unprime { get; set; } // O:    S_PRIME    Factor that denormalizes Kt', Kn', etc.
        public double Utime { get; set; } // T:      S_GEOM     Universal (Greenwich) standard time
        public double Zenetr { get; set; } /* T:     S_ZENETR   Solar zenith angle, no atmospheric correction (= ETR) */
        public double Zenref { get; set; } /* O:     S_REFRAC   Solar zenith angle, deg. from zenith, refracted */
    }



    public class SolPosTest
    {

        /*============================================================================
        *
        *    NAME:  stest00.c
        *
        *    PURPOSE:  Exercises the solar position algorithms in 'solpos.c'.
        *
        *        S_solpos
        *            INPUTS:     year, daynum, hour, minute, second, latitude,
        *                        longitude, timezone
        *
        *            OPTIONAL:   press   DEFAULT 1013.0 (standard pressure)
        *                        temp    DEFAULT   10.0 (standard temperature)
        *                        tilt    DEFAULT    0.0 (horizontal panel)
        *                        aspect  DEFAULT  180.0 (South-facing panel)
        *                        month   (if the S_DOY function is turned off)
        *                        day     ( "             "             "     )
        *
        *            OUTPUTS:    amass, ampress, azim, cosinc, coszen, day, daynum,
        *                        elevref, etr, etrn, etrtilt, month, prime,
        *                        sbcf, sretr, ssetr, unprime, zenref
        *
        *       S_init        (optional initialization for all input parameters in
        *                      the posdata struct)
        *           INPUTS:     struct posdata*
        *           OUTPUTS:    struct posdata*
        *
        *                     (Note: initializes the required S_solpos INPUTS above
        *                      to out-of-bounds conditions, forcing the user to
        *                      supply the parameters; initializes the OPTIONAL
        *                      S_solpos inputs above to nominal values.)
        *
        *      S_decode       (optional utility for decoding the S_solpos return code)
        *           INPUTS:     long int S_solpos return value, struct posdata*
        *           OUTPUTS:    text to stderr
        *
        *
        *        All variables are defined as members of the struct posdata
        *        in 'solpos00.h'.
        *
        *    Usage:
        *         In calling program, along with other 'includes', insert:
        *
        *              #include "solpos00.h"
        *
        *    Martin Rymes
        *    National Renewable Energy Laboratory
        *    25 March 1998
        *
        *    28 March 2001 REVISION:  SMW changed benchmark numbers to reflect the
        *                             February 2001 changes to solpos00.c
        *
        *----------------------------------------------------------------------------*/

        public string Main()
        {
            StringBuilder sb = new StringBuilder();

            Postdata pdat = new Postdata(); /* declare a posdata struct and a pointer for
	        //posdata pdat;
                                 it (if desired, the structure could be
                                 allocated dynamically with malloc) */
            ErrorCode retval; ; // to capture S_solpos return codes


            /**************  Begin demo program **************/

            //pdat = pd; // point to the structure for convenience

            /* Initialize structure to default values. (Optional only if ALL input
               parameters are initialized in the calling code, which they are not
               in this example.) */

            //S_init(pdat); /* C# Already done in constructor */

            /* I use Atlanta, GA for this example */

            pdat.Longitude = -84.43; // Note that latitude and longitude are
            pdat.Latitude = 33.65; //   in DECIMAL DEGREES, not Deg/Min/Sec
            pdat.Timezone = -5.0; /* Eastern time zone, even though longitude would
                                  suggest Central.  We use what they use.
                                  DO NOT ADJUST FOR DAYLIGHT SAVINGS TIME. */

            pdat.Year = 1999; // The year is 1999.
            pdat.Daynum = 203; /* July 22nd, the 203'rd day of the year (the
                                  algorithm will compensate for leap year, so
                                  you just count days). S_solpos can be
                                  configured to accept month-day dates; see
                                  examples below.) */

            /* The time of day (STANDARD time) is 9:45:37 */

            pdat.Hour = 9;
            pdat.Minute = 45;
            pdat.Second = 37;

            /* Let's assume that the temperature is 27 degrees C and that
               the pressure is 1006 millibars.  The temperature is used for the
               atmospheric refraction correction, and the pressure is used for the
               refraction correction and the pressure-corrected airmass. */

            pdat.Temp = 27.0;
            pdat.Press = 1006.0;

            /* Finally, we will assume that you have a flat surface facing southeast,
               tilted at latitude. */

            pdat.Tilt = pdat.Latitude; // Tilted at latitude
            pdat.Aspect = 135.0; // 135 deg. = SE

            sb.AppendLine("\n");
            sb.AppendLine("***** TEST S_solpos: *****\n");
            sb.AppendLine("\n");

            var solPos = new Solpos();

            retval = solPos.ExSolpos(ref pdat); // ExSolpos function call
            sb.AppendLine(solPos.S_decode(retval, ref pdat)); // ALWAYS look at the return code!

            /* Now look at the results and compare with NREL benchmark */

            sb.AppendLine("Note that your final decimal place values may vary\n");
            sb.AppendLine("based on your computer's floating-point storage and your\n");
            sb.AppendLine("compiler's mathematical algorithms.  If you agree with\n");
            sb.AppendLine("NREL's values for at least 5 significant digits, assume it works.\n\n");

            sb.AppendLine("Note that S_solpos has returned the day and month for the\n");
            sb.AppendLine("input daynum.  When configured to do so, S_solpos will reverse\n");
            sb.AppendLine("this input/output relationship, accepting month and day as\n");
            sb.AppendLine("input and returning the day-of-year in the daynum variable.\n");
            sb.AppendLine("\n");
            sb.AppendLine("NREL    -> 1999.07.22, daynum 203, retval 0, amass 1.335752,         ampress 1.326522\n");
            sb.AppendLine($"SOLTEST -> {pdat.Year:D}.{pdat.Month:D2}.{pdat.Day:D2}, daynum {pdat.Daynum:D}, retval {retval:D}, amass {pdat.Amass}, ampress {pdat.Ampress}\n");
            sb.AppendLine("NREL    -> azim 97.032875,        cosinc 0.912569,          elevref 48.409931\n");
            sb.AppendLine($"SOLTEST -> azim {pdat.Azim}, cosinc {pdat.Cosinc}, elevref {pdat.Elevref}\n");
            sb.AppendLine("NREL    -> etr 989.668518,       etrn 1323.239868,      etrtilt 1207.547363\n");
            sb.AppendLine($"SOLTEST -> etr {pdat.Etr}, etrn {pdat.Etrn}, etrtilt {pdat.Etrtilt}\n");
            sb.AppendLine("NREL    -> prime 1.037040,         sbcf 1.201910,         sunrise 347.173431\n");
            sb.AppendLine($"SOLTEST -> prime {pdat.Prime}, sbcf {pdat.Sbcf}, sunrise {pdat.Sretr}\n");
            sb.AppendLine("NREL    -> sunset 1181.111206,      unprime 0.964283,          zenref 41.590069\n");
            sb.AppendLine($"SOLTEST -> sunset {pdat.Ssetr}, unprime {pdat.Unprime}, zenref {pdat.Zenref}\n");



            /**********************************************************************/
            /* S_solpos configuration examples using the function parameter.

               Selecting a minimum of functions to meet your needs may result in
               faster execution.  A factor of two difference in execution speed
               exists between S_GEOM (the minimum configuration) and S_ALL (all
               variables calculated).  [S_DOY is actually the simplest and fastest
               configuration by far, but it only does the date conversions and bypasses
               all solar geometry.] If speed is not a consideration, use the default
               S_ALL configuration implemented by the call to S_init.

               The bitmasks are defined in S_solpos.h. */

            /* 1) Calculate the refraction corrected solar position variables */

            pdat.Function = Functions.S_REFRAC;


            /* 2) Calculate the shadow band correction factor */

            pdat.Function = Functions.S_SBCF;


            /* 3) Select both of the above functions (Note that the two bitmasks
                  are 'or-ed' together to produce the desired results): */

            pdat.Function = Functions.S_REFRAC | Functions.S_SBCF;


            /* 4) Modify the above configuration for accepting month and day rather
                  than day-of-year.  Note that S_DOY (which controls on the day-of-year
                  interpretation) must be inverted, then 'and-ed' with the other
                  function codes to turn the day-of-year OFF.  With the day-of-year
              bit off, S_solpos expects date input in the form of month and day. */

            pdat.Function = (Functions.S_REFRAC | Functions.S_SBCF) & ~Functions.S_DOY;
            pdat.Month = 7;
            pdat.Day = 22;

            /*    Also note that S_DOY is the only function that you should attempt
                  to clear in this manner: Other function bitmasks are a composite
                  of more than one mask, which represents an interdependency among
                  functions. Turning off unexpected bits will produce unexpected
                  results.  If in the course of your program you need fewer
                  parameters calculated, you should rebuild the function mask
                  from zero using only the required function bitmasks. */


            /* 5) Switch back to day-of-year in the above configuration by 'or-ing'
                  with the bitmask */

            pdat.Function |= Functions.S_DOY;
            pdat.Month = -99; // Now ignore ridiculous month and day
            pdat.Day = -99; // and overwrite them with correct values

            /*    ... and back to month-day again, etc.: */

            pdat.Function &= ~Functions.S_DOY;


            /* To see the effects of different configurations, move the above
               lines of code to just before the call to S_solpos and examine
               the results.  Note that some of the returned parameters are undefined
               if the configuration doesn't specify that they be calculated.  Thus,
               you must keep track of what you're calculating if you stray from the
               S_ALL default configuration. */


            /**********************************************************************/
            /* Looking at the S_solpos return code

               In the return code, each bit represents an error in the range of
               individual input parameters.  See the bit definition in S_solpos.h
               for the position of each error flag.

               To assure that none of your input variables are out of bounds, the
               calling program should always look at the S_solpos return code.  In
               this example, the function S_decode fulfills that mandate by examining
               the code and writing an interpretation to the standard error output.

               To see the effect of an out of bounds parameter, move the following
               line to just before the call to S_solpos: */

            pdat.Year = 99; // will S_solpos accept a two-digit year?

            /* This causes S_decode to output a descriptive line regarding the value
               of the input year. [This algorithm is valid only between years 1950 and
               2050; hence, explicit four-digit years are required. If your dates are
               in a two-digit format, S_solpos requires that you make a conversion
               to an explicit four-digit year.]

               S_decode (located in the solpos.c file) can serve as a template for
               building your own decoder to handle errors according to the
               requirements of your calling application. */


            /***********************************************************************/
            /* Accessing the individual functions */

            /* S_solpos was designed to calculate the output variables using the
               documented input variables.  However, as a matter of advanced
               programming convenience, the individual functions within S_solpos
               are accessible to the calling program through the use of the primative
               L_ masks (these are different from the composite S_ masks used
               above).  However, USE THESE WTTH CAUTION since the calling program
               must supply ALL parameters required by the function.  Because many of
               these variables are otherwise carefully created internally by
               S_solpos, the individual functions may not have bounds checking;
               hence your calling program must do all validation on the function
               input parameters. By the same reasoning, the return error code
               (retval) may not have considered all relevant input values, leaving
               the function vulnerable to computational errors or an abnormal end
               condition.

               As with the S_ masks above, the function variable is set to the
               L_ mask.  L_ masks may be ORed if desired.

               The solpos.h file contains a list of all output and transition
               variables, the reqired L_ mask, and all variables necessary for the
               calculation within individual functions.

               For example, the following code seeks only the amass value.  It calls
               only the airmass routine, which requires nothing but refracted zenith
               angle and pressure. Normally, the refracted zenith angle is a
               calculation that depends on many other functions within S_solpos.  But
               here using the L_ mask, we can simply set the refracted zenith angle
               externally and call the airmass function. */

            pdat.Function = Functions.L_AMASS | Functions.L_DOY; // call only the airmass function
            pdat.Press = 1013.0; // set your own pressure

            /* set up for the output of this example */
            sb.AppendLine("Raw airmass loop:\n");
            sb.AppendLine("NREL    -> 37.92  5.59  2.90  1.99  1.55  1.30  1.15  1.06  1.02  1.00\n");
            sb.Append("SOLTEST -> ");

            /* loop through a number of externally-set refracted zenith angles */
            for (pdat.Zenref = 90.0; pdat.Zenref >= 0.0; pdat.Zenref -= 10.0)
            {
                retval = solPos.ExSolpos(ref pdat); // call solpos                
                sb.Append($"{pdat.Amass,5:f2} "); // print out the airmass
                if (retval != 0) sb.AppendLine(solPos.S_decode(retval, ref pdat)); // retval may not be valid
            }
            sb.AppendLine("\n");

            return sb.ToString();
        }



    }


}
