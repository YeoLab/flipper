# NOTE TOO SELF: I should really save this for the time when I actually have the module set up so I do not need to walk people through all of this nonsense. 

# YEO-LAB internal example 
Hello. This is a short example for running Flipper as a member of the Yeo-lab partition on TSCC.

## Load up an interactive node.

```
srun -N 1 -c 1 -t 4:00:00 -p gold -q hcg-csd792 -A csd792 --mem 4G --pty /bin/bash
```

## Create a folder in scratch to save the output.
Run the command below (Replacing YOUR_USERNAME with your TSCC username) to create a folder to save the output of the 

```
mkdir /tscc/lustre/ddn/scratch/YOUR_USERNAME/Flipper_test
cd /tscc/lustre/ddn/scratch/YOUR_USERNAME/Flipper_test
```

## Load up the flipper module.

```
module load flipper
```

**IMPORTANT**: The only modules that should be loaded are the default modules (the ones that automatically load whenever you go on TSCC) and Skipper (Loading Skipper also automatically loads the singularity module). No other modules should be loaded, as this may confuse snakemake as to which environment to run the code in, causing an error. You can check what modules you have loaded using ```module list ```

## Copy the example config file from the module. 

Copy the config file from the skipper module folder to the scratch directory you just made. 

```
cp $/tscc/projects/ps-yeolab4/software/flipper/1.0.0/bin/flipper/example/yeo_lab_internal_example_config.yaml ./yeo_lab_internal_example_config.yaml
```

##  Adjust the yeo loab internal example config file. 
Adjust the WOKDIR input in your copy of `yeo_lab_internal_example_config.yaml` so that it points to the scratch directory you made in the first step. The config file can be adjsuted using any text editor you are comfortable with (vim, nano, jupyternotebooks, etc). 

No other changes to the config are necessary. 

## Run skipper. 

Now, simply replace YOUR_USERNAME in the command below again and run the following commands. 

```
unset SLURM_JOB_ID
   
snakemake -s $SKIPPER_HOME/bin/skipper/Skipper.py --configfile /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper2_test/yeo_lab_internal_example_config.yaml --profile $SKIPPER_HOME/bin/skipper/profiles/tscc2_snakemake9
```

**NOTE:** If problems occur in this initial run (or any run), please check the troubleshooting section below (especially point 2). If problems persist, please open up a github issue at [https://github.com/YeoLab/skipper](https://github.com/YeoLab/skipper)