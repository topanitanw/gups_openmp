#include <hpcc.h>
#include "RandomAccess.h"
#include <stdio.h>

/* Verification phase: local buckets to sort into */
#define BUCKET_SIZE 1024
#define SLOT_CNT 1
#define FIRST_SLOT 2

#define ZERO64B 0L

extern u64Int *Table;

void
HPCC_RandomAccessCheck(u64Int logTableSize,
					   u64Int logNumThreads,
					   u64Int TableSize,
					   u64Int NumUpdates,
					   s64Int *NumErrors)

{

  u64Int Ran;
  u64Int RanTmp;
  s64Int NextSlot;
  s64Int WhichPe;
  s64Int PeBucketBase;
  s64Int SendCnt;
  s64Int errors;
  int i;
  int j;
  s64Int PeCheckDone;
  int LocalAllDone =  HPCC_FALSE;
  int sAbort, rAbort;
  int numthreads, tid;
  u64Int *LocalBuckets;     /* buckets used in verification phase */
  u64Int *GlobalBuckets;    /* buckets used in verification phase */
  u64Int GlobalStartMyProc;
  u64Int LocalTableSize;

#pragma omp parallel private (tid, GlobalStartMyProc, PeCheckDone, PeBucketBase, Ran, WhichPe, NextSlot, j, RanTmp, i, GlobalBuckets)
 {

   numthreads = omp_get_num_threads();
   tid = omp_get_thread_num();

#pragma omp barrier

  GlobalBuckets = XMALLOC( u64Int, (numthreads*(BUCKET_SIZE+FIRST_SLOT)));
  sAbort = 0; if (! GlobalBuckets) sAbort = 1;

  if (sAbort > 0) {
    if (tid == 0) fprintf(stderr, "Failed to allocate memory for global buckets.\n");
    exit(1);
  }

  SendCnt = NumUpdates;

  LocalTableSize = TableSize / numthreads;

#pragma omp barrier

  GlobalStartMyProc = tid * LocalTableSize;

  Ran = starts (4 * GlobalStartMyProc);

  PeCheckDone = HPCC_FALSE;

#pragma omp barrier

  while(LocalAllDone == HPCC_FALSE){
    if (SendCnt > 0) {
      /* Initalize local buckets */

      PeBucketBase = tid * (BUCKET_SIZE+FIRST_SLOT);
      GlobalBuckets[PeBucketBase+SLOT_CNT] = FIRST_SLOT;
      GlobalBuckets[PeBucketBase+HPCC_DONE] = HPCC_FALSE;

      /* Fill local buckets until one is full or out of data */
      NextSlot = FIRST_SLOT;

      while(NextSlot != (BUCKET_SIZE+FIRST_SLOT) && SendCnt>0 ) {
        Ran = (Ran << 1) ^ ((s64Int) Ran < 0 ? POLY : 0);
        PeBucketBase = tid * (BUCKET_SIZE+FIRST_SLOT);
        NextSlot = GlobalBuckets[PeBucketBase+SLOT_CNT];
        GlobalBuckets[PeBucketBase+NextSlot] = Ran;
        GlobalBuckets[PeBucketBase+SLOT_CNT] = ++NextSlot;

	if (SendCnt != 0)
#pragma omp atomic
	  SendCnt--;
#pragma omp barrier
      }


    } /* End of sending loop */

    if (SendCnt == 0)
      GlobalBuckets[(tid*(BUCKET_SIZE+FIRST_SLOT))+HPCC_DONE] = HPCC_TRUE;

#pragma omp barrier

    LocalAllDone = HPCC_TRUE;

      if(PeCheckDone == HPCC_FALSE) {
        PeBucketBase = tid * (BUCKET_SIZE+FIRST_SLOT);
        PeCheckDone = GlobalBuckets[PeBucketBase+HPCC_DONE];
        for (j = FIRST_SLOT; j < GlobalBuckets[PeBucketBase+SLOT_CNT]; j ++) {
          RanTmp = GlobalBuckets[PeBucketBase+j];
#pragma omp critical
          Table[RanTmp & (TableSize-1)] ^= RanTmp;
        }
#pragma omp atomic
        LocalAllDone &= PeCheckDone;
      }

#pragma omp barrier

  }

#pragma omp master
  errors = 0;

  for (i=tid * LocalTableSize; i< (tid * LocalTableSize) + LocalTableSize; i++)

    {
    if (Table[i] != i + (LocalTableSize * tid) )
#pragma omp atomic
      errors++;
    }
 }

  *NumErrors = errors;

  free( GlobalBuckets );

  failed_globalbuckets:

  return;
}

