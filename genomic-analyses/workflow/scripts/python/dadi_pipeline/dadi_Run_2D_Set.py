'''
Usage: python dadi_Run_2D_Set.py
This is a modified version of the 'dadi_Run_Optimizations.py' script in which
we run optimizations for 2D comparisons for a large set of models that have been
made available as part of published works. These models are stored in the
Models_2D.py script, and will be called directly here. The user can delete or
comment out models to analyze a subset of the models available. 
This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_2D.py script, which
has all the model definitions.
General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.
 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688
Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 
Citations:
 If you use these scripts or the main diversification models for your work, please
 cite the following publication:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266
 
 If you use the additional diversification models or the island models set please cite 
 the following publication:
    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K.,
    Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. Sky, sea,
    and forest islands: diversification in the African leaf-folding frog
    Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae).
    Journal of Biogeography 45: 1781-1794. 
    doi: 10.1111/jbi.13365
        
 If you are interesting in contributing your models to this workflow, please email me!
-------------------------
Written for Python 2.7 and 3.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------
Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated September 2019
'''

import sys
import os
import numpy
import dadi
import pylab
import logging
import functools
from datetime import datetime
import Optimize_Functions
import Models_2D

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

# numpy.random.seed(42)
# dadi.Integration.timescale_factor = 1e-4

# Set up logger and send logs, stdout, and stderr to same file
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(snakemake.log[0], 'a')
fh_formatter = logging.Formatter('%(asctime)s %(levelname)s %(lineno)d:%(filename)s(%(process)d) - %(message)s')
fh.setFormatter(fh_formatter)
logger.addHandler(fh)
sys.stdout = open(snakemake.log[0], 'a')
sys.stderr = sys.stdout
print = functools.partial(print, flush = True)

#**************
snps = snakemake.input['sfs'][0]
#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["pop0", "pop1"]

#**************
#projection sizes, in ALLELES not individuals
proj = [30, 30]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#================================================================================
# Calling external 2D odels from the Models_2D.py script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:
 
 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, 
                        reps=None, maxiters=None, folds=None, in_params=None, 
                        in_upper=None, in_lower=None, param_labels=None, optimizer="log_fmin")
 
   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, 
          ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).
   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order
     optimizer: a string, to select the optimizer. Choices include: "log" (BFGS method), "log_lbfgsb" (L-BFGS-B method), 
                "log_fmin" (Nelder-Mead method), and "log_powell" (Powell's method).
Below, I give all the necessary information to call each model available in the
Models_2D.py script. I have set the optimization routine to be the same for each
model using the optional lists below, which are included as optional arguments for
each model. This particular configuration will run 4 rounds as follows:
Round1 - 10 replicates, maxiter = 3, fold = 3
Round2 - 20 replicates, maxiter = 5, fold = 2
Round3 - 30 replicates, maxiter = 10, fold = 2
Round4 - 40 replicates, maxiter = 15, fold = 1
If this script was run as is, each model would be called and optimized sequentially;
this could take a very long time. For your actual analyses, I strongly recommend
creating multiple scripts with only a few models each and running them
independently. It is also not a good idea to mix models from the Diversification Set
and the Island Set, as each was meant to be mutually exclusive.
'''


#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
prefix = '{0}{1}_{2}'.format(snakemake.params['prefix'], snakemake.wildcards.city, snakemake.wildcards.rep) 

#**************
#make sure to define your extrapolation grid size (based on your projections)
pts = [80,90,100]

#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

# Dictionary with parameters for all models
model_params_dict = {
        'no_div' : { 'params' : [''], 'lower' : [''], 'upper' : [''] },
        'no_div_bot' : { 'params' : ['nuB', 'nuF', 'TB', 'TF'], 'lower' : [0, 0, 0, 0], 'upper' : [20, 20, 20, 20] },
        'no_div_growth' : { 'params' : ['nu', 'T'], 'lower' : [0, 0], 'upper' : [20, 20] },
        'no_div_bot_growth' : { 'params' : ['nuB', 'nuG', 'TB', 'TG'], 'lower' : [0, 0, 0, 0], 'upper' : [20, 20, 20, 20] },
        'split_no_mig' : { 'params' : ['nu1', 'nu2', 'T'], 'lower' : [0, 0, 0], 'upper' : [20, 20, 20] },
        'split_sym_mig' : { 'params' : ['nu1', 'nu2', 'm', 'T'], 'lower' : [0, 0, 0, 0], 'upper' : [20, 20, 20, 20] },
        'split_asym_mig' : { 'params' : ['nu1', 'nu2', 'm12', 'm21', 'T'], 'lower' : [0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20] },
        'split_bot_urb_no_mig' : { 'params' : ['nu1', 'nu2', 'nu1B', 'nu1F', 'Tb', 'Ts'], 'lower': [0, 0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20, 20] },
        'split_bot_urb_sym_mig' : { 'params' : ['nu1', 'nu2', 'nu1B', 'nu1F', 'm', 'Tb', 'Ts'], 'lower' : [0, 0, 0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20, 20, 20] },
        'split_bot_urb_asym_mig' : { 'params' : ['nu1', 'nu2', 'nu1B', 'nu1F', 'm12', 'm21', 'Tb', 'Ts'], 'lower' : [0, 0, 0, 0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20, 20, 20, 20] },
        'split_bot_rur_no_mig' : { 'params' : ['nu1', 'nu2', 'nu2B', 'nu2F', 'Tb', 'Ts'], 'lower' : [0, 0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20, 20] },
        'split_bot_rur_sym_mig' : { 'params' : ['nu1', 'nu2', 'nu2B', 'nu2F', 'm', 'Tb', 'Ts'], 'lower' : [0, 0, 0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20, 20, 20] },
        'split_bot_rur_asym_mig' : { 'params' : ['nu1', 'nu2', 'nu2B', 'nu2F', 'm12', 'm21', 'Tb', 'Ts'], 'lower' : [0, 0, 0, 0, 0, 0, 0, 0], 'upper' : [20, 20, 20, 20, 20, 20, 20, 20] },
        'split_growth_urb_no_mig' : { 'params' : ['nu1', 's', 'T'], 'lower' : [0, 0, 0], 'upper' : [20, 0.99, 20] },
        'split_growth_urb_sym_mig' : { 'params' : ['nu1', 's', 'T', 'm'], 'lower' : [0, 0, 0, 0], 'upper' : [20, 0.99, 20, 20] },
        'split_growth_urb_asym_mig' : { 'params' : ['nu1', 's', 'T', 'm12', 'm21'], 'lower' : [0, 0, 0, 0, 0], 'upper' : [20, 0.99, 20, 20, 20] },
        'split_growth_rur_no_mig' : { 'params' : ['nu2', 's', 'T'], 'lower' : [0, 0, 0], 'upper' : [20, 0.99, 20] },
        'split_growth_rur_sym_mig' : { 'params' : ['nu2', 's', 'T', 'm'], 'lower' : [0, 0, 0, 0], 'upper': [20, 0.99, 20, 20] },
        'split_growth_rur_asym_mig' : { 'params' : ['nu2', 's', 'T', 'm12', 'm21'], 'lower' : [0, 0, 0, 0, 0], 'upper' : [20, 0.99, 20, 20, 20] }
        }

# Setup arguments for optimization routine
model = snakemake.wildcards.model
params = model_params_dict[model]['params']
lower = model_params_dict[model]['lower']
upper = model_params_dict[model]['upper']
param_labels = ','.join(params)
n_params = len(params)
model_func = getattr(Models_2D, model)

Optimize_Functions.Optimize_Routine(fs, pts, prefix, model, model_func, rounds, n_params, fs_folded=fs_folded,
                                    reps=reps, maxiters=maxiters, folds=folds, param_labels = param_labels,
                                    in_upper = upper, in_lower = lower)