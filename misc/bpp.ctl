* This is a BPP 4.0 control file !!!

          seed = -1

       seqfile = ../../neutral.25loci.4sp.txt
      Imapfile = ../../myImap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

  speciesdelimitation = 0     * fixed species tree

speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees; 2:user probs

  species&tree = 4  human  chimp  gorilla  orang
		 1	   1 	 	1      1
                 (((human, chimp), gorilla), orang);

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 25   * number of data sets in seqfile

     cleandata = 0      * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.008   * 2 500    # invgamma(a, b) for theta
      tauprior = 3 0.036   * 4 219 1  # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

      locusrate = 0       # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
       heredity = 0       # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)

      finetune = 1: .012 .003 .0001 .00005 .004 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 2000
      sampfreq = 10
       nsample = 10000

       *** Note: Make your window wider (144 columns) before running the program.
