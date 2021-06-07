/*-----------------------------------------------------------------------------

  RandomAccess Benchmark - UPC version - RandomAccess_UPC.c

  -- X1 ready version -- 10/28/04 --

  This benchmark is a UPC version of the RandomAccess code, developed by the
  High Performance Computing Laboratory at the George Washington University.

  Information on the UPC project at GWU is available at:

           http://upc.gwu.edu

  This benchmark is derived from the code from the High Performance Computing
  Challenge Benchmark Suite available at http://icl.cs.utk.edu/hpcc/

  UPC version:      F. Cantonnet    - GWU - HPCL (fcantonn@gwu.edu)
                    Y. Yao          - GWU - HPCL (yyy@gwu.edu)
                    T. El-Ghazawi   - GWU - HPCL (tarek@gwu.edu)

  Authors (HPCC):   D. Koester      - MITRE
                    B. Lucas        - USC / ISI
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program
    can be freely redistributed provided that you conspicuously and
    appropriately publish on each copy an appropriate referral notice to
    the authors and disclaimer of warranty; keep intact all the notices
    that refer to this License and to the absence of any warranty; and
    give any other recipients of the Program a copy of this License along
    with the Program.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA02111-1307 USA
-----------------------------------------------------------------------------*/

/* The benchmark is successfully completed if at least 99% of the table
 * is correctly built; it *is* OK, for example, to have a few errors
 * enter into the table due to rare race conditions between processors
 * updating the same memory location.
 *
 * It is expected that in a multi-processor version, one will start the
 * random number generator at N equally spaced spots along its cycle.
 * Each processor then steps the generator through its section of the
 * cycle.  (The starts routine provided can be used to start the
 * random number generator at any specified point.)
 *
 * The benchmark objective is to measure interprocessor bandwidth.  It
 * is ok to bucket sort locally and use a message passing protocol
 * between processors if this is faster than a global shared memory
 * approach.  More sophisticated optimizations which attempt to reorder
 * the loop so that all updates are local are not considered "in the
 * spirit" of the benchmark. */

#include <hpcc.h>
#include <omp.h>

/*#define LONG_IS_64BITS 1 */

#include "RandomAccess.h"
#include "string.h"

/* Number of updates to table (suggested: 4x number of table entries) */
#define NUPDATE (4 * TableSize)

/* Allocate main table (in global memory) */
u64Int *Table;
u64Int logTableSize, TableSize;
u64Int rank;
#define MAX_THREADS 16384
u64Int sh_temp[MAX_THREADS];

#define MAXJOBS 16384

void
RandomAccessUpdate(u64Int TableSize)
{
  s64Int i;
  static u64Int _ran[MAXJOBS];              /* Current random numbers */
  u64Int * const ran2 = _ran;
  int j;
  int tid;
  int numthreads;
  int numupdates;
  u64Int ran;
  double r1,rantime;
  u64Int LocalTableSize;
  omp_lock_t writelock;
  
  rantime = 0;

  /* Initialize main table */

  omp_init_lock(&writelock);

#pragma omp parallel private (tid, ran, i,r1, rantime)
 {
   tid = omp_get_thread_num();
   numthreads = omp_get_num_threads();
   numupdates = NUPDATE;



#pragma omp barrier

   LocalTableSize = TableSize / numthreads;

#pragma omp barrier

  for( i=tid * LocalTableSize  ; i< (tid * LocalTableSize) + LocalTableSize; i++ )
    Table[i] = i + (LocalTableSize * tid);

#pragma omp barrier
  /* Perform updates to main table.  The scalar equivalent is:
   *
   *     u64Int ran;
   *     ran = 1;
   *     for (i=0; i<NUPDATE; i++) {
   *       ran = (ran << 1) ^ (((s64Int) ran < 0) ? POLY : 0);
   *       table[ran & (TableSize-1)] ^= ran;
   *     }
   */

  ran = starts(4 * tid * LocalTableSize);
  rantime = 0;
#pragma omp barrier



#pragma omp for
  for(i = 0; i < numupdates; i++){
    

  ran = (ran << 1) ^ ((s64Int) ran < 0 ? POLY : 0);

  r1 = -RTSEC();
  omp_set_lock(&writelock);
  
  Table[(ran & (TableSize-1))] ^= ran;

  omp_unset_lock(&writelock);

  r1 += RTSEC();
  rantime += r1;

  }

  printf("%d critical and update time is %.6f\n", tid, rantime);

 }

 omp_destroy_lock(&writelock); 

}

int main(int argc, char **argv )
{
  double d_gups;
  int pow2_size;
  int r, i_failure;
  char *tmp;
  HPCC_Params p;

  if( argc == 2 )
    {
      sscanf( argv[1], "%d", &pow2_size );
      for( r=0, p.HPLMaxProcMem=sizeof(u64Int); 
	   r<pow2_size; r++, p.HPLMaxProcMem <<= 1 );
    }
  else
    p.HPLMaxProcMem = 20000000;

  tmp = p.outFname;
  strcpy(tmp, "hpccoutf.txt" );
  r = RandomAccess( &p, 1, &d_gups, &i_failure );

  return r;
}

int
RandomAccess(HPCC_Params *params, int doIO, double *GUPs, int *failure)
{
  s64Int i, j;
  u64Int temp;
  double cputime;               /* CPU time to update table */
  double realtime;              /* Real time to update table */
  double totalMem;
  FILE *outFile;
  int MYTHREAD;
  int numthreads;
  int tid;
  int numupdates;
  s64Int NumErrors;
  int logNumThreads;
  int PowerofTwo;
  u64Int localLogTableSize;
  u64Int localTableSize;
  

  outFile = fopen( params->outFname, "a" );

  if (! outFile)
    {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      exit(1);
    }

  /* calculate local memory per node for the update table */
  totalMem = params->HPLMaxProcMem;
  totalMem /= sizeof(u64Int);

  /* calculate the size of update array (must be a power of 2) */
  for (totalMem *= 0.5, logTableSize = 0, TableSize = 1;
       totalMem >= 1.0;
       totalMem *= 0.5, logTableSize++, TableSize <<= 1)
    ; /* EMPTY */


  localLogTableSize = logTableSize;
  localTableSize = TableSize;

#pragma omp parallel
  { 

#pragma omp master
      {
	numthreads =  omp_get_num_threads();
	fprintf( outFile, "Running on %d processors\n", numthreads );
	TableSize =  TableSize * numthreads;      

    for (i = 1, logNumThreads = 0; ; logNumThreads++, i <<= 1) {
      if (i == numthreads) {
	PowerofTwo = HPCC_TRUE;
	break;

      }
    
    }
    
    logTableSize = logTableSize + logNumThreads;

      }
  }

  fprintf( outFile, "Local table size   = 2^" FSTR64 " = " FSTR64 " words\n", 
	   (long long int)localLogTableSize, (long long int)localTableSize);

  fprintf( outFile, "Main table size   = 2^" FSTR64 " = " FSTR64 " words\n", 
	   (long long int)logTableSize, (long long int)TableSize);
  fprintf( outFile, "Number of updates = " FSTR64 "\n", (long long int)NUPDATE);

  Table = (u64Int *) XMALLOC (u64Int, TableSize);
  if (! Table)
    {
      if (doIO)
        {
          fprintf( outFile, "Failed to allocate memory for the update table (" FSTR64 ").\n", 
                            (long long int) TableSize);
          fclose( outFile );
        }
      exit(1);
    }

  /* Begin timing here */
  cputime = -CPUSEC();
  realtime = -RTSEC();

  RandomAccessUpdate( TableSize );

  /* End timed section */
  cputime += CPUSEC();
  realtime += RTSEC();

      /* make sure no division by zero */
      *GUPs = (realtime > 0.0 ? 1.0 / realtime : -1.0);
      *GUPs *= 1e-9*NUPDATE;

      /* Print timing results */

      fprintf( outFile, "CPU time used  = %.6f seconds\n", cputime);
      fprintf( outFile, "Real time used = %.6f seconds\n", realtime);
      fprintf( outFile, "%.9f Billion(10^9) Updates    per second [GUP/s]\n", *GUPs );


  /* Verification of results (in serial or "safe" mode; optional) */


      numupdates = NUPDATE;

      HPCC_RandomAccessCheck(logTableSize, logNumThreads, TableSize, numupdates, &NumErrors);
      
      fprintf( outFile, "Found " FSTR64 " errors in " FSTR64 " locations (%s).\n",
	       NumErrors, TableSize, (NumErrors <= 0.01*TableSize) ?
	       "passed" : "failed");

  return 0;
}

/* Utility routine to start random number generator at Nth step */
u64Int
starts(s64Int n)
{
  /* s64Int i, j; */
  int i, j;
  u64Int m2[64];
  u64Int temp, ran;

  while (n < 0)
    n += PERIOD;
  while (n > PERIOD)
    n -= PERIOD;
  if (n == 0)
    return 0x1;

  temp = 0x1;
  for (i=0; i<64; i++)
    {
      m2[i] = temp;
      temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
      temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
    }

  for (i=62; i>=0; i--)
    if ((n >> i) & 1)
      break;

  ran = 0x2;
  while (i > 0)
    {
      temp = 0;
      for (j=0; j<64; j++)
        if ((ran >> j) & 1)
          temp ^= m2[j];
      ran = temp;
      i -= 1;
      if ((n >> i) & 1)
        ran = (ran << 1) ^ ((s64Int) ran < 0 ? POLY : 0);
    }

  return ran;
}

