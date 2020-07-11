/*****************************************************************************
**                         CONV.C  version: 1.2                             **
******************************************************************************
** This program converts JPL-supplied unix binary files into a file readable**
** with C code under Linux and DOS on the Intel 80x86 based machines.       **
** It should also be readable with Fortran under both Linux and DOS.        **
******************************************************************************
**              This (1.2) version works with JPL ephemerides:              **
**              DE200, DE403, DE404, DE405 and DE406.                       **
******************************************************************************
**   This program reads unix binary file supplied by JPL ( DE405.UNX for    **
**   example) and creates binary file, in which all data from               **
**   the first file are copied WITH REVERSED ORDER OF BYTES                 **
**   in 4 byte integers and  8 byte double numbers.                         **
**   All other data and the overall structure ( record length,              **
**   variable order etc.) remains the same.                                 **
**   It is also a possibility of creating the cinverted binary file only    **
**   for a given time subinterval of the whole ephemeris.                   **
******************************************************************************
**    YOU SHOULD MANUALLY ADJUST THIS SOURCE FILE BEFORE COMPILED AND RUN   **
******************************************************************************
** Originally written by: Piotr A. Dybczynski, e-mail:dybol@amu.edu.pl      **
******************************************************************************
** Created: July  2, 1997 by PAD ** Last modified: July 23, 1997 by PAD.    **
*****************************************************************************/

#include<stdio.h>
#include<math.h>

void reverse( unsigned char t[], int n);

/****************************************************************************/
/* UNCOMMENT ONE AND ONLY ONE OF THE FOLLOWING DENUM DEFINITIONS: */

/*#define DENUM 200*/
/*#define DENUM 403*/
/*#define DENUM 404*/
#define DENUM 405
/*#define DENUM 406*/
/****************************************************************************/

#if   DENUM==200
#define KSIZE 1652
#elif DENUM==403 || DENUM==405
#define KSIZE 2036
#elif DENUM==404 || DENUM==406
#define KSIZE 1456
#endif

#define NRECL 4
#define RECSIZE (NRECL*KSIZE)
#define NCOEFF (KSIZE/2)

/* we use long int instead of int for DOS-Linux compatibility: in either case
   there are 4 bytes per integer variable !   */

/* Declarations in comments comes from ephemeris reading software.
   They are replaced with char arrays for byte order reversing.    */

struct rec1{
         char ttl[3][84];
         char cnam[400][6];
         /*double ss[3];*/
         unsigned char ss[24];
         /*long int ncon;*/
         unsigned char ncon[4];
         /*double au;*/
         unsigned char au[8];
         /*double emrat;*/
         unsigned char emrat[8];
         /*
         long int ipt[12][3];
         long int numde;
         long int lpt[3];
         */
         unsigned char other[160];
       };
 struct{
         struct rec1 r1;
         char spare[RECSIZE-sizeof(struct rec1)];
       } R1;

 struct rec2{
         /*double cval[400];*/
         unsigned char cval[3200];
       };
 struct{
         struct rec2 r2;
         char spare[RECSIZE-sizeof(struct rec2)];
       } R2;

union {
        double ss[3];
        unsigned char bb[24];
      } B;


FILE *F1,*F2;

void main(void)
{
  int i;
  unsigned char buff[10000];
  double t1,t2;

  puts("CONV: file conversion program for the JPL export ephemerides\n");
  printf("Enter original JPL binary file name: ");
  scanf(" %s",buff);

  F1=fopen(buff,"rb");
  if(F1==NULL){
                printf("Cannot open %s file for reading, aborted!\n",buff);
                return;
              }

  printf("Output file name: ");
  scanf(" %s",buff);

  F2=fopen(buff,"wb");
  if(F2==NULL){
                printf("Cannot open %s file for writting, aborted!\n",buff);
                fclose(F1);
                return;
              }


/* reading first two records (with constants, epochs, pointers etc. */

  fread(&R1,RECSIZE,1,F1);
  fread(&R2,RECSIZE,1,F1);

  for(i=0;i<3;++i) reverse(R1.r1.ss+8*i,8);
  for(i=0;i<24;++i) B.bb[i]=R1.r1.ss[i];
  reverse(R1.r1.ncon,4);
  reverse(R1.r1.au,8);
  reverse(R1.r1.emrat,8);
  for(i=0;i<40;++i) reverse(R1.r1.other+4*i,4);

  for(i=0;i<400;++i) reverse(R2.r2.cval+8*i,8);

  /* Decide what interval is desired....              */
  printf("Enter Start Epoch (in JED, zero if all file): ");
  scanf(" %lf",&t1);
  if(t1!=0.0)
    {
      if(t1>B.ss[0])
        {
         t1 -= 2440400.5;
         t1 = B.ss[2]*floor(t1/B.ss[2])+2440400.5;
         B.ss[0]=t1;
        }
      printf("Enter Final Epoch (in JED): ");
      scanf(" %lf",&t2);
      if(t2<B.ss[1])
        {
         t2 -= 2440400.5;
         t2 = B.ss[2]*ceil(t2/B.ss[2])+2440400.5;
         B.ss[1]=t2;
        }
      for(i=0;i<16;++i) R1.r1.ss[i]=B.bb[i];
    }
  else
    {
      t1=0.0;
      t2=999999999.0;
    }
/* writing first two records converted and adjusted for new interval */

  fwrite(&R1,RECSIZE,1,F2);
  fwrite(&R2,RECSIZE,1,F2);

/* main loop of reading, converting and writting ephemeris binary data */

  while( !feof(F1) )
       {
         i=fread(buff,RECSIZE,1,F1);
         if(i!=1) break;
         reverse(buff,8);
         reverse(buff+8,8);
         for(i=0;i<16;++i) B.bb[i]=buff[i];
         if(B.ss[0] >= t2) break;
         if(B.ss[1] <= t1) continue;
         for(i=2;i<NCOEFF;++i) reverse(buff+i*8,8);
         fwrite(buff,RECSIZE,1,F2);
       }

  fclose(F1);
  fclose(F2);
}
/**************************************************************************/

/* subroutine reverse reverses the order of n nbytes starting from t */

void reverse( unsigned char t[], int n)
{
  unsigned char buff[10];
  int i;

  for(i=0;i<n;++i) buff[n-i-1] = t[i];
  for(i=0;i<n;++i) t[i] = buff[i];
}