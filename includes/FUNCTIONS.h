/* 
    Author: Marina Reggiani Guzzo (The University of Manchester)

    This file contains all the functions common to my codes, so 
    I can always refer to them throughout the codes where they 
    might be necessary. Which guarantees that I am using the same 
    function among my codes. So if I have to fix the definition of
    any function, the updated one will automatically be used in 
    all codes.
*/

#ifndef FUNCTIONS
#define FUNCTIONS

// ROOT input
#include "TChain.h"
#include "vector"
#include "TH1.h"
#include "TVector3.h"

// checkout input
#include "eval.h"
#include "kine.h"
#include "pfdump.h"
#include "pot.h"
#include "tagger.h"
#include "pfeval.h"

// my libraries
#include "myFunctions.C"       // functions to get input data

#endif