1. The [01_species_trees](https://github.com/suz11001/Tripartition/tree/main/sim_data/01_species_trees) directory contains a 100 different species tree.   
2. The [02_gene_trees](https://github.com/suz11001/Tripartition/tree/main/sim_data/02_gene_trees) directory contains a single gene tree for each of the 100 species tree. Note that the pruned gene tree must be placed in a seperate directory as shown.  
3. The [03_sub-gene_trees](https://github.com/suz11001/Tripartition/tree/main/sim_data/03_sub-gene_trees) contains three sub-directories, low, medium, and high. These directories correspond to sub-gene transfer rate of 0.2, 0.4, and 0.6, respectively. Note that the pruned sub-gene tree and leaf mapping must be placed in seperate directories as showns. 

The commands used to generate these data sets are explained [here](https://github.com/suz11001/Tripartition/blob/main/sim_data/SimulatedData_CommandsUsed.pdf)

4. The [04_sequences](https://github.com/suz11001/Tripartition/tree/main/sim_data/04_sequences) directory contains sequences 6 sub-directories. Each sub-directory represents the percent of the sub-gene length i.e. `subgene_len_20perc` means all the sub-genes are 20% of the total gene length. Within each sub-gene length directory, there are nine sub-directories, they represent the following parameters:  
   a. 2000bps  = gene length is 2000 bps   
   b. 500bps   = gene length is 500 bps  
   c. control  = default setting i.e. gene length = 1000 bps, sub-gene length = 40% (400 bps), sub-gene transfer rate = 0.4, substitution rate = 0.5  
   d. domain_ht  = sub-gene transfer rate is 0.6  
   f. domain_lt  = sub-gene transfer rate is 0.4  
   g. scaled_0.1 = substitution rate is 0.1  
   h. scaled_1  = substitution rate is 1  
   i. scaled_2  = substitution rate is 2  
   j. scaled_5 = substitution rate is 5  

   Note that only one parameter is varied at a time i.e. the default values are used except for the one parameter that we are varying.

The command used to generate the sequences is following (for the baseline dataset):   
`python3 partial.py -z $random -dl 400 -gl 600 -ap -s 0.5 /path/to/pruned/subgene/tree /path/to/pruned/gene/tree /path/to/mapping/file/`
   
