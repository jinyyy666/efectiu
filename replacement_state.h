#ifndef REPL_STATE_H
#define REPL_STATE_H
 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This file is distributed as part of the Cache Replacement Championship     //
// workshop held in conjunction with ISCA'2010.                               //
//                                                                            //
//                                                                            //
// Everyone is granted permission to copy, modify, and/or re-distribute       //
// this software.                                                             //
//                                                                            //
// Please contact Aamer Jaleel <ajaleel@gmail.com> should you have any        //
// questions                                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cassert>
#include "utils.h"
#include "crc_cache_defs.h"
#include <iostream>

using namespace std;

// define the associavity for the dead block sampler
#define SAMPLER_ASSOC 12

// define the sampling gap for the LLC
#define SAMPLE_GAP 64

// LLC shared cache size: 4MB and there are 4096 sets
// define total number of sample sets for the sampler
#define TOTAL_SAMPLE_SETS 4096/SAMPLE_GAP

// define total bits of partial tags and PC
#define PARTIAL_BITS 15

// define the total number of prediction tables (skew optimization)
#define TOTAL_PREDICT_TABLES 3

// define the replacement policy for the sampler, 0 -- LRU 1 -- random
#define SAMPLER_REPLACEMENT_POLICY 0

// Replacement Policies Supported
typedef enum 
{
    CRC_REPL_LRU        = 0,
    CRC_REPL_RANDOM     = 1,
    CRC_REPL_CONTESTANT = 2
} ReplacemntPolicy;

// Replacement State Per Cache Line
typedef struct
{
    UINT32  LRUstackposition;
    bool dead; // One bit metadata to indicate whether the block is dead;
    // CONTESTANTS: Add extra state per cache line here

} LINE_REPLACEMENT_STATE;

// Cache block for sampler
struct m_block{
    UINT32 lru_position;
    Addr_t partial_PC;  // 15-bit partial PC
    Addr_t partial_tag; // 15-bit partial tag
    bool prediction;    // one-bit information for dead-block prediction for sampler itself
    bool valid;
    Addr_t PC;          // full PC

    m_block (void){
        partial_PC = 0;
	partial_tag = 0;
	PC = 0;
        valid = false;
	prediction = false;
    }
};

// Set for sampler:
struct m_set{
    m_block blocks[ SAMPLER_ASSOC ];
    bool set_valid;
   
    m_set (void){
        set_valid = false;
	for(UINT32 i=0; i < SAMPLER_ASSOC; ++i)
	    blocks[i].lru_position = i;
    }
};

// Sampler for the dead block prediction
struct SAMPLER{
    m_set * sets;
    UINT32 nsets;
    UINT32 assoc;
    INT32 replacement_policy;
  
    SAMPLER( void ){}
};  


// Predictor (2-bit saturated counter table)
typedef struct{
    UINT32 cnt;
} PREDICTOR;


// The implementation for the cache replacement policy
class CACHE_REPLACEMENT_STATE
{
public:
    LINE_REPLACEMENT_STATE   **repl;
    SAMPLER *sampler;
    PREDICTOR **predictor;
  private:

    UINT32 numsets;
    UINT32 assoc;
    UINT32 replPolicy;

    COUNTER mytimer;  // tracks # of references to the cache

    // CONTESTANTS:  Add extra state for cache here
    bool dead_block_detection; // indication of finding any dead blocks in LLC
    UINT32 blreplPolicy;       // baseline replacement policy for LLC under db sampling 
    Addr_t prev_PC;            // the previous PC which don't corresponding to prefetch
  public:
    ostream & PrintStats(ostream &out);

    // The constructor CAN NOT be changed
    CACHE_REPLACEMENT_STATE( UINT32 _sets, UINT32 _assoc, UINT32 _pol );

    INT32 GetVictimInSet( UINT32 tid, UINT32 setIndex, const LINE_STATE *vicSet, UINT32 assoc, Addr_t PC, Addr_t paddr, UINT32 accessType );
    void   UpdateReplacementState( UINT32 setIndex, INT32 updateWayID, const LINE_STATE *currLine, 
                                   UINT32 tid, Addr_t PC, UINT32 accessType, bool cacheHit );

    void   SetReplacementPolicy( UINT32 _pol ) { replPolicy = _pol; } 
    void   IncrementTimer() { mytimer++; } 


    ~CACHE_REPLACEMENT_STATE(void);

  private:
    
    void   InitReplacementState();
    INT32  Get_Random_Victim( UINT32 setIndex );
    INT32  Get_My_Victim( 
			UINT32 setIndex, const LINE_STATE *vicSet, 
			UINT32 assoc, Addr_t PC, Addr_t paddr, UINT32 accessType ); 
    INT32  Get_LRU_Victim( UINT32 setIndex );
    void   UpdateLRU( UINT32 setIndex, INT32 updateWayID );
    void   UpdateMyPolicy( UINT32 setIndex, INT32 updateWayID, const LINE_STATE *currLine,  
			    Addr_t PC, UINT32 accessType, bool cacheHit);
    void   AccessSampler(UINT32 sampleIndex, Addr_t tag, Addr_t PC, UINT32 accessType);
    void   MoveToMRU(m_set & samplerSet, UINT32 way);
    void   UpdatePredictionTable(Addr_t PC, bool increase);
    bool   GetPrediction(Addr_t PC, const INT32 const_thresh = 8);
    void   PlaceSamplerSet(m_set & samplerSet, UINT32 i, Addr_t tag, Addr_t PC);
    UINT32 hash(Addr_t key, Addr_t a_prime, Addr_t b_prime, Addr_t p_prime);
};

#endif
