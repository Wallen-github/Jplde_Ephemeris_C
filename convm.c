/*****************************************************************************
**                        CONVM.C  version: 1.2                             **
******************************************************************************
** This program converts MULTIPLE JPL-supplied unix binary files into       **
** a single binary file readable with C code under Linux and DOS on the     **
** Intel 80x86 based machines. It should also be readable with Fortran      **
** under both Linux and DOS. Tests for various fortran compilers in progress**
******************************************************************************
**              This (1.2) version works with JPL ephemerides:              **
**              DE200, DE403, DE404, DE405 and DE406.                       **
******************************************************************************
**   This program reads a series of standard unix binary file supplied      **
**   by JPL ( UNXP1600.405, UNXP1650.405, ... for EXAMPLE ) and creates     **
**   a single, binary ephemeris  file, in which all data from               **
**   the original files were copied WITH REVERSED ORDER OF BYTES            **
**   in 4 byte integers and  8 byte double numbers.                         **
**   All other data and the overall structure ( record length,              **
**   variable order etc.) remains the same.                                 **
**   It is also a possibility of creating the converted binary file only    **
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
#define DSTART 1600
#define DSTOP  2150
#define DSTEP  50
#define STEP   32.
#elif DENUM==404 || DENUM==406
#define KSIZE 1456
#define DSTART -3000
#define DSTOP   2700
#define DSTEP  300
#define STEP   64.
#endif

#define NRECL 4
#define RECSIZE (NRECL*KSIZE)
#define NCOEFF (KSIZE/2)

/* we use long int instead of int for DOS-Linux compatibility: in either case
   there are 4 bytes per integer variable !   */

/* Declarations in comments come from ephemeris reading software.
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
  int i,first,firstblock,fnum;
  unsigned char buff[10000];
  double t1,t2,te1,tlast,tend;

  puts("CONVM: multiple file conversion program for the JPL export ephemerides\n");
   /* Decide what interval is desired....              */
  printf("Enter Start Epoch (in JED, zero if whole ephemeris): ");
  scanf(" %lf",&t1);
  if(t1!=0.0)
    {
      t1 -= 2440400.5;
      t1 = STEP*floor(t1/STEP)+2440400.5;
      printf("Enter Final Epoch (in JED): ");
      scanf(" %lf",&t2);
      t2 -= 2440400.5;
      t2 = STEP*ceil(t2/STEP)+2440400.5;
    }
  else
    {
      t1=0.0;
      t2=999999999.0;
    }

  printf("\nOutput file name: ");
  scanf(" %s",buff);

  F2=fopen(buff,"wb");
  if(F2==NULL){
                printf("Cannot open %s file for writting, aborted!\n",buff);
                fclose(F1);
                return;
              }

  first=1;
  firstblock=1;

  for(fnum=DSTART;fnum<=DSTOP;fnum+=DSTEP)
     {

/************ YOU SHOULD ADJUST THE FOLLOWING LINE TO YOUR ENVIRONMENT*******/

       sprintf(buff,"D:\\de405\\unx%c%04d.%03d",(fnum<0?'m':'p'),(fnum<0?-1*fnum:fnum),DENUM);

/****************************************************************************/

       F1=fopen(buff,"rb");
       if(F1==NULL)
         {
           printf("File: %s is not available.\n",buff);
           continue;
         }

/**** reading first two records from each file (constants, epochs, etc.) ****/

       fread(&R1,RECSIZE,1,F1);
       fread(&R2,RECSIZE,1,F1);

/**** converting epochs to compare with selected interval *******************/
       for(i=0;i<3;++i) reverse(R1.r1.ss+8*i,8);
       for(i=0;i<16;++i) B.bb[i]=R1.r1.ss[i];

/* write once the first two records converted. They are identical in all
   binary files except two variables: starting and final epochs.
   The first record will be overwritten at the end of conversion with those
   epochs adjusted to the ephemeris file created */

       if(first)
         {
           reverse(R1.r1.ncon,4);
           reverse(R1.r1.au,8);
           reverse(R1.r1.emrat,8);
           for(i=0;i<40;++i) reverse(R1.r1.other+4*i,4);

           for(i=0;i<400;++i) reverse(R2.r2.cval+8*i,8);

           fwrite(&R1,RECSIZE,1,F2);
           fwrite(&R2,RECSIZE,1,F2);
           first=0;
         }

/**** checking whether this file contens useful data ************************/

       if(B.ss[0] >= t2) break;      /* all data for the subinterval copied */
       if(B.ss[1] <= t1) continue;      /* this file is for previous epochs */

/* main loop of reading, converting and writting ephemeris binary data */
       while( fread(buff,RECSIZE,1,F1)==1)
            {
/*** fist two variables stored in each block (8 bytes each) are starting
     and final epochs for that block of chebyshev polinomial coefficients ***/

              reverse(buff,8);
              reverse(buff+8,8);
              for(i=0;i<16;++i) B.bb[i]=buff[i];

              if(B.ss[0] >= t2) break;
              if(B.ss[1] <= t1) continue;
              if(firstblock)
                {
                  te1=B.ss[0];                /* saving first epoch written */
                  firstblock=0;
                  tlast=te1;
                  tend=B.ss[1];
                }
              else
                {
                  if(B.ss[0]==tlast) continue;  /* allow for double blocks */
                  if(B.ss[0]!=tend)
                    {
                      puts("\nAborted becouse of the gap in binary data.\n");
                      fnum=9999;
                      break;
                    }
                }

              for(i=2;i<NCOEFF;++i) reverse(buff+i*8,8);
              fwrite(buff,RECSIZE,1,F2);
              tlast=B.ss[0];
              tend=B.ss[1];
            }
       fclose(F1);
     }

/*** all binary data were converted and stored in the new file **************/
/*** now overwrite the first record with corrected epochs *******************/

  fflush(F2);
  fseek(F2,0,SEEK_SET);                /* point to the begining of the file */
  B.ss[0]=te1;                                       /* restore first epoch */
  B.ss[1]=tend;                                       /* restore last epoch */
  for(i=0;i<16;++i) R1.r1.ss[i] = B.bb[i];   /* copy them into first record */
  /* reverse the rest of values */
  reverse(R1.r1.ncon,4);
  reverse(R1.r1.au,8);
  reverse(R1.r1.emrat,8);
  for(i=0;i<40;++i) reverse(R1.r1.other+4*i,4);

  fwrite(&R1,RECSIZE,1,F2);     /* ovewrite first record with actual epochs */

  fclose(F2);
}
/****************************************************************************/

/* subroutine reverse reverses the order of n nbytes starting from t */

void reverse( unsigned char t[], int n)
{
  unsigned char buff[10];
  int i;

  for(i=0;i<n;++i) buff[n-i-1] = t[i];
  for(i=0;i<n;++i) t[i] = buff[i];
}