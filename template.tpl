//Parameters for the coalescence simulation program : simcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NPOP0
NPOP1
NPOP2
//Samples sizes and samples age 
$sample_number_of_location1
$sample_number_of_location2
0
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG_01 0
MIG_10 0 0
MIG_20 MIG_21 0
//Migration matrix 1
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
6 historical event
TDIV 0 2 1 RESIZE0 0 1
TDIV 1 2 1 RESIZE1 0 1
TDIV 2 2 0 RESIZE2 0 1
TDIV 1 1 0 1 0 1
TDIV 0 0 0 1 0 1
TDIV 2 2 0 1 0 1  
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0.4773642 2e-8
