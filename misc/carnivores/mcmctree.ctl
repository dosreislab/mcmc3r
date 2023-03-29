          seed = -1
       seqfile = morph_mol.aln
      treefile = carnivores.tree
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 2
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = 'B(37.3, 66.0, 0.025, 0.025)'
       TipDate = 1 1  * TipDate (1) & time unit

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5  * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0 0.001 * birth, death, sampling
   alpha_gamma = 1 1         * gamma prior for alpha

   rgene_gamma = 2 10  * gamma prior for overall rates for genes
  sigma2_gamma = 2 2   * gamma prior for sigma^2 (for clock=2 or 3)

         print = 2
	burnin = 500
      sampfreq = 1 
       nsample = 8000
