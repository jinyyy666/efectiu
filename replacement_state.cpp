#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>
#include <map>
#include <iostream>

using namespace std;

#include "replacement_state.h"


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

/*
** This file implements the cache replacement state. Users can enhance the code
** below to develop their cache replacement ideas.
**
*/

static Addr_t primes[TOTAL_PREDICT_TABLES][3] = 
{
    {997,359,57},
    {521,131,73},
    {337,677,17}
}; // a table of primes for hashing function


////////////////////////////////////////////////////////////////////////////////
// The replacement state constructor:                                         //
// Inputs: number of sets, associativity, and replacement policy to use       //
// Outputs: None                                                              //
//                                                                            //
// DO NOT CHANGE THE CONSTRUCTOR PROTOTYPE                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
CACHE_REPLACEMENT_STATE::CACHE_REPLACEMENT_STATE( UINT32 _sets, UINT32 _assoc, UINT32 _pol )
{

    numsets    = _sets;
    assoc      = _assoc;
    replPolicy = _pol;

    mytimer    = 0;

    InitReplacementState();
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// The function prints the statistics for the cache                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
ostream & CACHE_REPLACEMENT_STATE::PrintStats(ostream &out)
{

    out<<"=========================================================="<<endl;
    out<<"=========== Replacement Policy Statistics ================"<<endl;
    out<<"=========================================================="<<endl;

    // CONTESTANTS:  Insert your statistics printing here
    
    return out;

}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function initializes the replacement policy hardware by creating      //
// storage for the replacement state on a per-line/per-cache basis.           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void CACHE_REPLACEMENT_STATE::InitReplacementState()
{
    // Create the state for sets, then create the state for the ways

    repl  = new LINE_REPLACEMENT_STATE* [ numsets ];


    // ensure that we were able to create replacement state:

    assert(repl);

    // Create the state for the sets
    for(UINT32 setIndex=0; setIndex<numsets; setIndex++) 
    {
        repl[ setIndex ]  = new LINE_REPLACEMENT_STATE[ assoc ];

        for(UINT32 way=0; way<assoc; way++) 
        {
            // initialize stack position (for true LRU)
            repl[ setIndex ][ way ].LRUstackposition = way;
	    // initalize dead indication bit (for dead block sampler)
	    repl[ setIndex ][ way ].dead = false;
        }
    }

    if (replPolicy != CRC_REPL_CONTESTANT) return;

    // Set the baseline replacement policy for LLC, 0 -- LRU, 1 -- random
    blreplPolicy = 0;

    // Set the previous pc:
    prev_PC = 0;

    // Create the sampler, which is basically a cache

    sampler = new SAMPLER();

    // ensure that we were able to create replacement state:

    assert(sampler);

    // Create the sets for the samplers:
    sampler->sets = new m_set[TOTAL_SAMPLE_SETS];
    sampler->nsets = TOTAL_SAMPLE_SETS;
    sampler->assoc = SAMPLER_ASSOC;
    sampler->replacement_policy = SAMPLER_REPLACEMENT_POLICY;


    // Create the predictor, which is basically a 2-bit saturated counter table:
    
    predictor = new PREDICTOR* [TOTAL_PREDICT_TABLES];
    
    assert(predictor);
    
    // Initialize the counters in the table as zero:
    for(UINT32 table = 0; table < TOTAL_PREDICT_TABLES; ++table){
        predictor[ table ] = new PREDICTOR[ numsets ];
	for(UINT32 row = 0; row < numsets; ++row)
	    predictor[ table ][ row ].cnt = 0;
    }

    
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is called by the cache on every cache miss. The input        //
// argument is the set index. The return value is the physical way            //
// index for the line being replaced.                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
INT32 CACHE_REPLACEMENT_STATE::GetVictimInSet( UINT32 tid, UINT32 setIndex, const LINE_STATE *vicSet, UINT32 assoc, Addr_t PC, Addr_t paddr, UINT32 accessType ) {
    // If no invalid lines, then replace based on replacement policy
    if( replPolicy == CRC_REPL_LRU ) 
    {
        return Get_LRU_Victim( setIndex );
    }
    else if( replPolicy == CRC_REPL_RANDOM )
    {
        return Get_Random_Victim( setIndex );
    }
    else if( replPolicy == CRC_REPL_CONTESTANT )
    {
        // Contestants:  ADD YOUR VICTIM SELECTION FUNCTION HERE
        return Get_My_Victim (setIndex, vicSet, assoc, PC, paddr, accessType);
    }

    // We should never here here

    assert(0);
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is called by the cache after every cache hit/miss            //
// The arguments are: the set index, the physical way of the cache,           //
// the pointer to the physical line (should contestants need access           //
// to information of the line filled or hit upon), the thread id              //
// of the request, the PC of the request, the accesstype, and finall          //
// whether the line was a cachehit or not (cacheHit=true implies hit)         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void CACHE_REPLACEMENT_STATE::UpdateReplacementState( 
    UINT32 setIndex, INT32 updateWayID, const LINE_STATE *currLine, 
    UINT32 tid, Addr_t PC, UINT32 accessType, bool cacheHit )
{

    // What replacement policy?
    if( replPolicy == CRC_REPL_LRU ) 
    {
        UpdateLRU( setIndex, updateWayID );
    }
    else if( replPolicy == CRC_REPL_RANDOM )
    {
        // Random replacement requires no replacement state update
    }
    else if( replPolicy == CRC_REPL_CONTESTANT )
    {
        // Contestants:  ADD YOUR UPDATE REPLACEMENT STATE FUNCTION HERE
        // Update my policy according to the paper:
        UpdateMyPolicy( setIndex, updateWayID, currLine, PC, accessType, cacheHit );
         
    }
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//////// HELPER FUNCTIONS FOR REPLACEMENT UPDATE AND VICTIM SELECTION //////////
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function finds the LRU victim in the cache set by returning the       //
// cache block at the bottom of the LRU stack. Top of LRU stack is '0'        //
// while bottom of LRU stack is 'assoc-1'                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
INT32 CACHE_REPLACEMENT_STATE::Get_LRU_Victim( UINT32 setIndex )
{
	// Get pointer to replacement state of current set

	LINE_REPLACEMENT_STATE *replSet = repl[ setIndex ];
	INT32   lruWay   = 0;

	// Search for victim whose stack position is assoc-1

	for(UINT32 way=0; way<assoc; way++) {
		if (replSet[way].LRUstackposition == (assoc-1)) {
			lruWay = way;
			break;
		}
	}

	// return lru way

	return lruWay;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function finds a random victim in the cache set                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
INT32 CACHE_REPLACEMENT_STATE::Get_Random_Victim( UINT32 setIndex )
{
    INT32 way = (rand() % assoc);
    
    return way;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function implements the LRU update routine for the traditional        //
// LRU replacement policy. The arguments to the function are the physical     //
// way and set index.                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void CACHE_REPLACEMENT_STATE::UpdateLRU( UINT32 setIndex, INT32 updateWayID )
{
	// Determine current LRU stack position
	UINT32 currLRUstackposition = repl[ setIndex ][ updateWayID ].LRUstackposition;

	// Update the stack position of all lines before the current line
	// Update implies incremeting their stack positions by one

	for(UINT32 way=0; way<assoc; way++) {
		if( repl[setIndex][way].LRUstackposition < currLRUstackposition ) {
			repl[setIndex][way].LRUstackposition++;
		}
	}

	// Set the LRU stack position of new line to be zero
	repl[ setIndex ][ updateWayID ].LRUstackposition = 0;
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function implements the dead block sampling replacement policy        //
// The arguments to the function are set index, set tag, associtivity,        //
// PC, physical memory address and the accessType                             //
// way and set index.                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

INT32 CACHE_REPLACEMENT_STATE::Get_My_Victim( 
        UINT32 setIndex, const LINE_STATE *vicSet, 
        UINT32 assoc, Addr_t PC, Addr_t paddr, UINT32 accessType ) 
{
        // Implement the bypassing policy here, if predicted dead, simply ignore!
        // Try to use the different threshold for bypassing!
        // Also ignore the writeback because of its irregular activity:
        if(GetPrediction(PC, 8) || accessType == ACCESS_WRITEBACK)
	    return -1;

        // Get pointer to the replacement state of the current set 
        // to find out whether or not there is block in LLC is predicted to be dead
        LINE_REPLACEMENT_STATE *replSet = repl[ setIndex ]; 

	// Find out if there is any dead blocks in LLC:
	UINT32 llc_way;
	for(llc_way = 0; llc_way < assoc; ++llc_way){
	    // always pickup the first predicted dead block
	    if(replSet[llc_way].dead)
	        break;
	}
	
	// assume no dead block is detected:
	dead_block_detection = false;
	
	// no block is predicted dead, use lru as default for LLC
	if(llc_way == assoc){
	    return Get_LRU_Victim( setIndex );
	}
	else{
            // if dead block is detected, return its way
	    dead_block_detection = true;
	    return (INT32)llc_way;
	}
}

void CACHE_REPLACEMENT_STATE::UpdateMyPolicy( 
    UINT32 setIndex, INT32 updateWayID, const LINE_STATE *currLine, 
    Addr_t PC, UINT32 accessType, bool cacheHit) 
{
  
         // ignore the writeback because of its irregular activity:
         if(accessType == ACCESS_WRITEBACK)
             return;      

	 if(accessType == ACCESS_PREFETCH){
	     // previous PC is needed bc all prefetch share
	     // the same PC!
	     //prev_PC = PC;
	 }

         // Get pointer to replacement state of the current set
         // in order to update the death indication bit in the LLC
         LINE_REPLACEMENT_STATE *replSet = repl[ setIndex ];
         
	 // Update the death indication bit for LLC according to the prediction table
	 // Use the 15-bit partial PC to index the prediction table
	 if(cacheHit && accessType != ACCESS_PREFETCH)
	     replSet[updateWayID].dead = GetPrediction(PC);
	 else if(!cacheHit)
	     replSet[updateWayID].dead = false;

	 // remember to update the lru stack position for LLC under LRU policy
	 // there are two cases: 
	 // 1) hit at LLC
	 // 2) misses at LLC and the victim block is selected by LRU
	 if(blreplPolicy == 0){
	     // case 1) if hit, update the lru stack position given the updateWayID
	     if(cacheHit)
	         UpdateLRU( setIndex, updateWayID );
	     else{
	         // case 2) if no dead block is detected in LLC during miss
	         // then update LRU
	         //if(!dead_block_detection){
		     UpdateLRU( setIndex, updateWayID );
		     //}  
	     }
	 }    
	 
	 // Need to update its prediction bit of the sampler set if accessed
	 // Also need to update the prediction table 
         // Note that I ignore the prefetch because all prefetches have the same PC!
	 if(setIndex % SAMPLE_GAP == 0){
 	     AccessSampler(setIndex/SAMPLE_GAP,currLine[updateWayID].tag,PC,accessType);
	 }

}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is to access the sampler and determine hit or miss           //
// If hit, decrease the prediction table entry index by partial PC            //
// this PC can be the incoming PC or PC stored in the sampler block           //
// If capacity/conflict miss, increase the prediction table entry in a         //
// similar way.                                                               //
// The arguments to the functions are sampler set index, tag of the           //
// incoming block and PC made this access                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void CACHE_REPLACEMENT_STATE::AccessSampler(UINT32 sampleIndex, Addr_t tag, Addr_t PC, UINT32 accessType)
{
    // Get the pointer to the set in the sampler:
    assert(sampleIndex < sampler-> nsets && sampleIndex >= 0);
    m_set & samplerSet = sampler->sets[sampleIndex];
    INT32 s_assoc = sampler->assoc;

    // mask to extract the 15-bit tag coming from the LLC and PC made this access:
    Addr_t mask = (1<<PARTIAL_BITS) - 1;
    Addr_t fullPC = PC;
    tag &= mask;
    PC &= mask;

    // Access the set to find out whether or not the tag is matched
    INT32 i = 0;
    for(i = 0; i < s_assoc; ++i){
        // if hit, get the prediction for the block and decrease the prediction table
        if(samplerSet.blocks[i].valid && samplerSet.blocks[i].partial_tag == tag){
	    // Use the previous PC to do the update: 
	    // if the access type is prefetch, please skip the update of the tab!
	    if(accessType != ACCESS_PREFETCH){
	        UpdatePredictionTable(samplerSet.blocks[i].partial_PC, false);

		samplerSet.blocks[i].prediction = GetPrediction(PC);
		samplerSet.blocks[i].partial_PC = PC;
		samplerSet.blocks[i].PC = fullPC;
	    }

            // please remember to update its LRU position as well if needed!
	    if(sampler->replacement_policy == 0)
	        MoveToMRU(samplerSet, i);

	    // NOTE: you can choose PC/ PC of the last instruct record in this block
	    //UpdatePredictionTable(PC, false);

	    return;
	}
    }
    
    // if miss
    bool set_valid = samplerSet.set_valid;
    if(!samplerSet.set_valid){
        // is there any invalid block left ?
        for(i = 0; i < s_assoc; ++i)
	    if(samplerSet.blocks[i].valid == false) break;
	
	if(i == s_assoc){
	    // mark the sampler set as having only valid blocks so 
	    // we don't search it again!
	    samplerSet.set_valid = true;
	    set_valid = true;
	}
	// at this point i indicates an invalid block
	// or assoc if there is no invalid block
    }

    // find out the dead block in the sampler, if none is found, i should be assoc:
    if(set_valid){
        INT32 j;
        for(j = 0; j < s_assoc; ++j){
	    // always find the first dead block
	  if(samplerSet.blocks[j].prediction){
	        break;
	  }
	}
	i = j;
    }
	
    if(i == s_assoc){
        // this is the case when no dead block is found
        assert(set_valid);
	
	// lru as the default replacement policy for sampler
	if(sampler->replacement_policy == 0){
	    i = s_assoc - 1;
	    assert(i>= 0 && i < s_assoc);

	    // use the evicted block signature (partial PC) to increase the table
	    UpdatePredictionTable(samplerSet.blocks[i].partial_PC, true);
	    

	    MoveToMRU(samplerSet, i);

	    // place the partial tag and partial PC into the sampler
	    // set the valid bit as true and prediction 
	    //if(accessType != ACCESS_PREFETCH)
	        PlaceSamplerSet(samplerSet, 0, tag, PC);

	    // if you place a faked PC, indicate that the block does not need to be updated
		samplerSet.blocks[0].PC = fullPC;
	}
	else{ // random policy as the default replacement policy for sampler
	    i = rand() % s_assoc;
	    
	    // use the evicted block signature (parital PC) to increase the table
	    UpdatePredictionTable(samplerSet.blocks[i].partial_PC, true);	    

	    PlaceSamplerSet(samplerSet, i, tag, PC);

	}
    }
    else{
        // this is the case 1) when dead block is found or 2) the block is invalid!
        assert((samplerSet.blocks[i].prediction == true && samplerSet.blocks[i].valid == true) 
	       || (samplerSet.blocks[i].prediction == false && samplerSet.blocks[i].valid == false));

	// need to update the prediction table if the block is found dead!
	if(samplerSet.blocks[i].valid){
	    UpdatePredictionTable(samplerSet.blocks[i].partial_PC, true);
	}    
	
	//MoveToMRU(samplerSet, i);

	PlaceSamplerSet(samplerSet, i, tag, PC);
	    
	samplerSet.blocks[0].PC = fullPC;
    }

}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is to move a block of the sampler to the MRU position        //
// The arguments to the function are sampler set and                          //
// index of the block to be moved                                             //
// samplerSet will be changed after this function                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void CACHE_REPLACEMENT_STATE::MoveToMRU(m_set & samplerSet, UINT32 way)
{
    assert(way >= 0 && way < SAMPLER_ASSOC);
    
    INT32 j; 
    m_block b = samplerSet.blocks[way];

    for(j = way; j >= 1; --j) samplerSet.blocks[j] = samplerSet.blocks[j-1];

    samplerSet.blocks[0] = b;

}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is to update the prediction table indexed by the hashed PC   //
// The arguments to the function are partial PC and                           //
// a flag to indicate whether increase/decrease the table                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void CACHE_REPLACEMENT_STATE::UpdatePredictionTable(Addr_t PC, bool increase)
{

}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is to generate a prediction based on the prediction table    //
// The argument to the function is the partial PC                             //
// The function returns the prediction results true->dead, false->not dead    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

bool CACHE_REPLACEMENT_STATE::GetPrediction(Addr_t PC, const INT32 const_thresh)
{
    INT32 sum = 0;
    for(UINT32 table = 0; table < TOTAL_PREDICT_TABLES; ++table){
        UINT32 entry = hash(PC, primes[table][0], primes[table][1], primes[table][2]);
	assert(entry >=0 && entry < numsets);
	assert(predictor[table][entry].cnt >=0 && predictor[table][entry].cnt < 4);
	
	sum += predictor[table][entry].cnt;
    }
    return sum >=  const_thresh ? true : false;
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is to place a block into the sampler                         //
// The arguments to the function are sampler set and                          //
// index of the updated way, partial tag, and partial PC                      //
// samplerSet will be changed after this function                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void CACHE_REPLACEMENT_STATE::PlaceSamplerSet(m_set & samplerSet, UINT32 i, Addr_t tag, Addr_t PC)
{
    assert(i >= 0 && i < SAMPLER_ASSOC);
    
    // set the valid bit as true and prediction as false 
    assert( tag < 1<<PARTIAL_BITS && tag >= 0 && PC < 1<<PARTIAL_BITS && PC >= 0);
    samplerSet.blocks[i].partial_tag = tag;
    samplerSet.blocks[i].partial_PC = PC;
    samplerSet.blocks[i].valid = true;
    samplerSet.blocks[i].prediction = GetPrediction(PC);

}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This function is to use three prime number to do a hashing                 //
// The arguments to the function are the partial PC as the key                //
// and three prime numbers                                                    //
// The function returns a table entry (0 ~ numsets)                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

inline
UINT32 CACHE_REPLACEMENT_STATE::hash( 
  Addr_t key, Addr_t a_prime, Addr_t b_prime, Addr_t p_prime)
{

}


CACHE_REPLACEMENT_STATE::~CACHE_REPLACEMENT_STATE (void) {
}
