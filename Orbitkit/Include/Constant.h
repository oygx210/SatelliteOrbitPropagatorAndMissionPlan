#if !defined(AFX_DEFINES_H_INCLUDED_)
#define AFX_DEFINES_H_INCLUDED_

#include <math.h>

// math functions                                                             
#define MIN(A, B)		((A) < (B) ? (A) : (B))
#define MAX(A, B)		((A) > (B) ? (A) : (B))
#define SQR(x)			((x)*(x))
#define fixangle(a)		((a) - 360.0 * (floor((a) / 360.0)))  
#define torad(d)		((d) * (PI / 180.0))                     
#define todeg(d)		((d) * (180.0 / PI))                     
#define RADIAN(d)		((d) * (PI / 180.0))                     
#define DEGREE(d)		((d) * (180.0 / PI))                     

/******************************************************************************/
/*                                                                            */
/* general numerical constants                                                */
/*                                                                            */
/******************************************************************************/

//注：长度单位均为:km
#ifndef RS
	#define RS 696000
#endif

#ifndef RE
	#define RE 6378.137
#endif

#ifndef RM
	#define RM 1738.2
#endif

#ifndef DIS_SE
	#define  DIS_SE 1.496e+8
#endif

#ifndef DIS_ME
	#define  DIS_ME 3.844e+5
#endif

#ifndef GM
	#define  GM 398600.5
#endif

#ifndef PI
	#define PI 3.1415926535897932384626433832795
#endif

#ifndef M_PI
	#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef TWOPI
	#define TWOPI 6.283185307179586476925286766559
#endif

#ifndef J2
	#define J2 0.00108263
#endif

#ifndef J3
	#define J3 -0.00000254
#endif

#ifndef J4
	#define J4 -0.00000161
#endif

#define FOURPI                (4.0*PI)

#ifndef HALFPI
#define HALFPI                (0.5*PI)
#endif

#define THREEHALFPI           (1.5*PI)

#define ZERO                  0.0
#define ONEPPB                1.0e-9
#define ONEPPM                1.0e-6
#define TWOPPM                2.0e-6
#define ONETHIRD              (1.0/3.0)
#define TWOTHIRDS             (2.0/3.0)
#define THREEHALFS            (3.0/2.0)
#define ONE                   1.0
#define ONEMEG                1.0e6
#define TWOMEG                2.0e6

/******************************************************************************/
/*                                                                            */
/* numerical constants for unit conversions                                   */
/*                                                                            */
/******************************************************************************/
#define CRH                   (24.0/TWOPI)         /* convert rad into hours  */
#define CRS                   (86400.0/TWOPI)      /* convert rad into sec    */
#define CRD                   (360.0/TWOPI)        /* convert rad into deg    */
#define CRAM                  (21600.0/TWOPI)      /* convert rad into arcmin */
#define CRAS                  (1296000.0/TWOPI)    /* convert rad into arcsec */
#define CDR                   (TWOPI/360.0)        /* convert deg into rad    */
#define CAMR                  (TWOPI/21600.0)      /* convert arcmin into rad */
#define CASR                  (TWOPI/1296000.0)    /* convert arcsec into rad */
#define CRREV                 (1.0/TWOPI)          /* convert rad into rev    */

#define HALFDEG               (0.5*CDR)            /* half a deg [rad]        */

#define MPD                   1440.0               /* minutes per day         */
#define MPD2                  (MPD*MPD)            /* (minutes per day)^2     */
#define MPD3                  (MPD2*MPD)           /* (minutes per day)^3     */
#define SPD                   86400.0              /* seconds per day         */

#define CKMM                  1000.0               /* convert km to m         */
#define CMKM                  1.0e-3               /* convert m to km         */
#define CKMNM                 0.539956804          /* convert km to naut. mil.*/
#define CHZKHZ                1.0e-3               /* convert Hz to kHz       */
#define CKHZMHZ               1.0e-3               /* convert kHz to MHz      */
#define CHZMHZ                1.0e-6               /* convert Hz to MHz       */
#define CKHZHZ                1.0e+3               /* convert kHz to Hz       */
#define CMHZHZ                1.0e+6               /* convert MHz to Hz       */

/******************************************************************************/
/*                                                                            */
/* numerical constants describing the Earth's orbit and figure                */
/*                                                                            */
/* EARTHRADIUS and EARTHFLAT are from: "Astronomical Almanac", 1991, p. K13   */
/*                                                                            */
/******************************************************************************/

#define EARTHSMA              149597892.0          /* 1 AU [km]               */
#define EARTHRADIUS           6378.137             /* equatorial radius [km]  */
#define EARTHECCEN            0.01675104           /* Earth's orbit eccentr.  */
#define EARTHFLAT             (1.0/298.257222)     /* geoid model parameters  */

/******************************************************************************/
/*                                                                            */
/* numerical constants describing the apparent size of the Sun and the        */
/* proximity limit for calculating transits across the solar disk             */
/*                                                                            */
/******************************************************************************/

#define SUNRADIUS             695980.0             /* equatorial radius [km]  */
#define SUNDISKRAD            (16.0*CAMR)          /* Sun disk radius [rad]   */
#define SUNPROX               48.0                 /* Sun prox limit [arcmin] */

/******************************************************************************/
/*                                                                            */
/* numerical constants describing the motions within the solar system         */
/*                                                                            */
/******************************************************************************/

#define JULCENT               36525.0              /* mean solar days / jcy   */
#define TROPCENT              36524.219879         /* mean solar days / cy    */
#define TROPYEAR              (TROPCENT/100.0)     /* mean solar days / year  */
#define JULDAT1900            2415020.0            /* Julian date of 1900.0   */
#define JULDAT1950            2433282.423          /* Julian date of 1950.0   */
#define JULDAT2000            2451545.0            /* Julian date of 2000.0   */
#define SIDSOLAR              1.002737909350       /* sidereal rotation rate  */
#define SIDRATE               (TWOPI/SPD*SIDSOLAR) /* [rad/s]                 */

/******************************************************************************/
/*                                                                            */
/* physical constants                                                         */
/*                                                                            */
/* GM (and KEPLER) are from: "Astronomical Almanac", 1991, p. K13             */
/*                                                                            */
/******************************************************************************/

#define CVAC                  2.99792458e5         /* speed of light [km/s]   */

#define GM                    398600.5		       /* [km^3/s^2]              */
#define GMSGP                 398600.7995          /* value used in SGP model */

#define KEPLER                42241.09773          /* GM^(1/3)*(SPD/2PI)^(2/3)*/
#define KEPLERSGP             42241.10831          /* value used in SGP model */
                                                   /* with a [km] and T [min] */
/* Julian date at standard epoch */
#define  J2000             2451545.0        
#define  TO_CENTURIES	   36525.0
#define  SEC_IN_DAY        86400.0
#define  HRS_IN_RADIAN     3.819718634205
#define  DEG_IN_RADIAN     57.2957795130823
/* flattening of earth, 1/298.257 */
#define  FLATTEN           0.003352813   
/* equatorial radius of earth, meters */
#define  EQUAT_RAD         6378137.0    
/* 1 AU in meters */
#define  ASTRO_UNIT        1.4959787066e11 

/*  Astronomical constants  */

#define epoch       2444238.5      /* 1980 January 0.0 */

/*  Constants defining the Sun's apparent orbit  */

#define elonge      278.833540     /* Ecliptic longitude of the Sun
                                      at epoch 1980.0 */
#define elongp      282.596403     /* Ecliptic longitude of the Sun at
                                      perigee */
#define eccent      0.016718       /* Eccentricity of Earth's orbit */
#define sunsmax     1.495985e8     /* Semi-major axis of Earth's orbit, km */
#define sunangsiz   0.533128       /* Sun's angular size, degrees, at
                                      semi-major axis distance */

/*  Elements of the Moon's orbit, epoch 1980.0  */

#define mmlong      64.975464      /* Moon's mean longitude at the epoch */
#define mmlongp     349.383063     /* Mean longitude of the perigee at the
                                      epoch */
#define mlnode      151.950429     /* Mean longitude of the node at the
                                      epoch */
#define minc        5.145396       /* Inclination of the Moon's orbit */
#define mecc        0.054900       /* Eccentricity of the Moon's orbit */
#define mangsiz     0.5181         /* Moon's angular size at distance a
                                      from Earth */
#define msmax       384401.0       /* Semi-major axis of Moon's orbit in km */
#define mparallax   0.9507         /* Parallax at distance a from Earth */
#define synmonth    29.53058868    /* Synodic month (new Moon to new Moon) */
#define lunatbase   2423436.0      /* Base date for E. W. Brown's numbered
                                      series of lunations (1923 January 16) */

/*  Properties of the Earth  */

#define earthrad    6378.137        /* Radius of Earth in kilometres */
#define EarthFlat (1/298.257222)            /* Earth Flattening Coeff. */

#ifndef OMEGAE
	#define OMEGAE 7.292116E-05
#endif

#define WM_USER_CREATEORBIT WM_USER+101
#define WM_USER_DELETEORBIT WM_USER+102
#define WM_USER_CREATEORBITSS WM_USER+103
#define WM_USER_DELETEORBITSS WM_USER+104
#define WM_USER_SHOWCONTROLBAR WM_USER+105
#define WM_USER_RELOADMODEL WM_USER+106

#define NEWCREATEDSAT 0
#define SGP4SDP4SAT   1

#endif