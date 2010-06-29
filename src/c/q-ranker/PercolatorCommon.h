/*******************************************************************************
 * Q-ranker unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 * and Marina Spivak (spivak.marina@gmail.com) at
 * NEC Labs of America.
 *******************************************************************************/
#ifndef PERCOLATOR_COMMON_H_
#define PERCOLATOR_COMMON_H_
#ifndef PERCOLATOR_C_INTERFACE_H_
#ifdef __cplusplus
extern "C" {
#endif

/** Number of target and decoy sets that external program will hand over to percolator.
 * The value should correspond to the number of sequence databases that have been searched.
 * Percolators validation strategy will be the same as for the stand alone version given the
 * corresponding number of sqt-files as input. */
typedef enum {TWO_SETS=2,THREE_SETS,FOUR_SETS} NSet;

/** Identifying which set the PSM belongs to*/
typedef enum {TARGET=0,DECOY1,DECOY2,DECOY3} SetType;
#ifdef __cplusplus
}
#endif
#endif
#endif
