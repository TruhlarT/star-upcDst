#ifndef RUNDEF_H
#define RUNDEF_H

#include <vector>

const bool runMAINANA = false; // true
const bool runEMBEDDING = false; // central embedding study
const bool runVERTEXSTUDY = true; // take quite long
const bool runTOFQA = false; // true Designed to be run on single run i.e. 1 job for 1 run
const bool runTRIGEFF = false; //true
const bool runFULLZB = false; // true
const bool runELASTICANA = false;
const bool runRPMCANA = false;
// ALIGNMENT is designed to be run separetly i.e. only ALIGNMENT 
// also is is designed to be run on single run i.e. 1 job for 1 run
const bool runALIGNMENT = false; 

const bool runStudy[] = { runMAINANA, runVERTEXSTUDY, runEMBEDDING, runTOFQA, runTRIGEFF, runFULLZB, 
            runELASTICANA, runRPMCANA, runALIGNMENT}; // must be in the same order as in enum in Util 

#endif //ifndef RUNDEF_H
