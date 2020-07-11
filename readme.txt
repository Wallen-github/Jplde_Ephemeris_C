*****************************************************************************
*****      C version software for the JPL planetary ephemerides.        *****
*****************************************************************************

version 1.2, June 23, 1997 (this README.txt file was updated April 19, 2007)

Piotr A. Dybczynski, Astronomical Observatory of the A. Mickiewicz Universty,
Sloneczna 36, 60-286 Poznan, POLAND.
e-mail: dybol@amu.edu.pl

*****************************************************************************
File list:
   Size   date     time   file  
  10218 2007-04-20 10:38 asc2eph.c
   6405 2007-04-20 10:38 conv.c
   9098 2007-04-20 10:38 convm.c
   1892 2007-04-20 10:39 jplbin.h
   5759 2007-04-20 10:47 readme.txt
  33421 2007-04-20 10:39 testeph.c
 409808 1998-06-03 00:00 testpo1.405
****************************************************************************

CHANGES FROM VERSION 1.1: added support for DE405 and DE406,
                          ASC2EPH.C can work with incoplete source file set,
                          added CONVM.C - new program for converting and
                                          merging multiple binary files,
                          bug in libration calculations fixed,
                          some minor changes and bug fixes

******************************************************************************
 The original JPL datafiles, FORTRAN code and documentation can be found at:
              http://ssd.jpl.nasa.gov/?planet_eph_export
******************************************************************************

This version (1.2) works with DE200, DE403, DE404, DE405 and DE406.

This file describes public domain software for using and manipulating
JPL export planetary ephemerides, written in C language.

In original JPL export packages you can find several FORTRAN source files
for reading and testing ephemerides.

Here you can find C versions of ASC2EPH and TESTEPH programs and an additional
two programs: CONV and CONVM, for converting original JPL suplied unix binary
files to DOS/LINUX binaries, without using any additional source files. Please
note, that we assume you do not rename any of the original, JPL supplied files.

ASC2EPH.C and TESTEPH.C are manual translations of the original JPL
FORTRAN code, slightly changed and adapted for C language.
They are designed to keep FORTRAN-based structure of the binary ephemeris
files so they can be shared between C and Fortran programmers working both
under DOS and LINUX on PC-based machines. ASC2EPH uses original JPL supplied
source (ASCII) ephemeris files. There is no need for editing them to
replace 'D' with 'E' in exponent parts of double precision numbers.

*****************************************************************************
IT IS NECESSARY TO ADJUST MANUALLY SOURCE FILES before you compile and run
them. First, look at the file: JPLBIN.H, used by ASC2EPH.c and TESTEPH.C
It contains the definition of DENUM, JPL ephemeris number. Some variables
obtain their values depending on DENUM (eg. RECSIZE - record size ).
Second, look at lines with fopen() calls in ASC2EPH.C, TESTEPH.C and CONV.C
Adjust them for your environment, adding path before filenames when necessary.
******************************************************************************

Program CONV.C is an additional tool for obtaining DOS/LINUX, FORTRAN/C
binary files with JPL ephemerides DIRECTLY from JPL supplied binary unix
files (*.unx). It simply reverses the order of bytes in each integer (4 byte)
and double (8 byte) number stored in the single JPL unix binary file found on
the original CD-ROM for example.
Additionally it can produce the subset of the whole ephemeris - you have to
specify starting and final JED epochs.

New program (v.1.2): CONVM.C is designed for simultaneous merging and
converting original JPL supplied binary UNIX files (you can obtain them from
an anonymous FTP site: ftp://navigator.jpl.nasa.gov ) into a single, binary
ephemeris file, readable with FORTRAN and C under both DOS and Linux. It works
similarly to CONV.C. You can specify smaller interval of time than the original
ephemeris time span.

All source files are rich of comments and ( I hope ) therefore are
self explanatory.

On July 8, 1997 JPL realeased corrected file TESTPO.405, without fictious 
"test points" for 5-th and 6-th components of nutation. You can find this 
file at: ftp://navigator.jpl.nasa.gov/ephem/export/test_data. There is a new 
file: testpo.406 there!

PLEASE NOTE that due to including of nutation and librations in the test points
file (GREAT!) you can expect differences greater than 1.0e-13 when using 
TESTEPH program. Some values of librations have more than three digits before 
decimal point and we do not have 17 significant digits on PC-based machines!
JPL have changed (June 23 version of TESTEPH.F) the formula for the difference
calculation from absolute (del = |xc-xo|) to relative: del = |xc-xo|/(|xc|+|xo|).
I have applied this correction in my software and found that now almost all del 
values are greater than 1e-13, generating software alert. So I decided to 
comment out this modification. Instead I used the old formula and I found that 
for all points in TESTPO.405 I have obtained del < 1e-13. This happened even for 
libration value for JD=2524897.5: 19440.0251465376386. My calculated value was: 
19440.0251465376386 ! I have 18 correct signifficant digits ! 
I have not switched to long double type, I still use standard 8 byte double 
precision arithmetic ;-).

All sugestions, bug reports or questions may be directed to the author:

Piotr A. Dybczynski, Internet e-mail: dybol@amu.edu.pl

Please send me a word if you use my software with success also.

April 19, 2007.

