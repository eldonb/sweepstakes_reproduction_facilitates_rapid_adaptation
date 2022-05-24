/* 2N diploid individuals in the population at the start; so 4N allele copies */
const double GLOBAL_CONST_N = 500000. ;
const double GLOBAL_CONST_A = 4. * GLOBAL_CONST_N ;
const double GLOBAL_CONST_I = 2. * GLOBAL_CONST_N ;
const unsigned int GLOBAL_CONST_II = static_cast<unsigned int>( GLOBAL_CONST_I) ;
const unsigned int GLOBAL_CONST_CUTOFF_ONE = 1000000 ;
const unsigned int GLOBAL_CONST_CUTOFF_TWO = 1000000 ;
const double GLOBAL_CONST_PSI_ONE = static_cast<double>(GLOBAL_CONST_CUTOFF_ONE) ;
const double GLOBAL_CONST_PSI_TWO = static_cast<double>(GLOBAL_CONST_CUTOFF_TWO) ;
const double GLOBAL_CONST_ALPHA_ONE = 0.75 ;
const double GLOBAL_CONST_ALPHA_TWO = 3.0 ;
const double GLOBAL_CONST_SELECTION = 0.5 ;
const unsigned int GLOBAL_CONST_BOTTLENECK = 10000 ;
const double GLOBAL_CONST_PROBABILITY_BOTTLENECK = 0.1 ;
const double GLOBAL_CONST_EPSILON =  0.1 ;
/* 1./static_cast<double>(GLOBAL_CONST_N) ; */
const int GLOBAL_CONST_NUMBER_EXPERIMENTS = 200000 ;
/* DRA : 1e4 experiments */
/* DRB : 1e4 experiments */
/* DRC : 1e6 experiments */
/* b  1e2 ; p 0.1 ; e 0 ; s 0.5 */
/* trajectoryDRC1366735.txt ;  trajectoryDRC9362053.txt */
/* DRD : 1e6 experiments */
/* b 1e4 ; p 0.1 ; e 0 ; s 0.5 */
/* 6 trajectoryDRD2970107.txt
     2 trajectoryDRD3395395.txt
     4 trajectoryDRD6857234.txt
     6 trajectoryDRD743686.txt
    10 trajectoryDRD9826111.txt
    28 total
*/
/*  DRE : 1e6 experiments
  b 1e2 ; p 0.01 ; e 0.1 ; s 0.5 
    2 trajectoryDRE1643468.txt
    4 trajectoryDRE3450340.txt
    1 trajectoryDRE6195849.txt
    4 trajectoryDRE8042068.txt
   11 total
*/
/*
  DRF : 1e4 experiments
  b 1e4 ; p 0.01 ; e 0.1 ; s 0.5 
1 trajectoryDRF7985480.tx
 */
/*
1e6 experiments
 b 1e2 ; p 0.1 ; e 0.1 ; s 0.5 
  1 trajectoryDRG1871443.txt
   2 trajectoryDRG5912666.txt
   3 total
 */
/*  1e6 experiments
  b 1e4 ; p 0.1 ; e 0.1 ; s 0.5 
  3 trajectoryDRH2423565.txt
    1 trajectoryDRH3015242.txt
    2 trajectoryDRH3678393.txt
    2 trajectoryDRH4574213.txt
    1 trajectoryDRH847514.txt
    9 total
 */
