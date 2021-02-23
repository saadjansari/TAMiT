#!/usr/bin/env python

import os
import sys
import pdb
import matlab.engine
import argparse
import re
import pdb
import subprocess
import numpy as np

'''
Name : RunSimulSingleCell.py
Description : Runs multiple fits on Summit by distirbuting them across multiple cores and nodes as needed.
Usage : RunSimulSingleCell.py summit
'''

def parse_args():
    # parse comand line arguments

    parser = argparse.ArgumentParser( prog='SimulSingleCell.py', description='Launcher for multiple matlab single cell fits on summit') 
    
    # Fit option: Run fits
    parser.add_argument('-F', '--fit', action='store_true',  
            help='launch fit jobs via sbatch')

    # Analyze flag: Run analysis 
    parser.add_argument('-A', '--analyze', action='store_true', 
            help='Analyze the cells either post fit or in specified directory if fit flag false.')

    # Analyze directory: directory for analysis, used when fit flag is not enabled
    parser.add_argument('--dir', type=str, 
            help='Directory for analysis, used when the fit flag is disabled.')

    # DryRun option: Run normal and create jobs but don't launch the jobs
    parser.add_argument('-n', '--dryrun', action='store_true',  
            help='run normal except dont launch the jobs')

    opts = parser.parse_args()

    # ensure dir is specified if certain conditions are met
    if not opts.fit and opts.analyze and not opts.dir:
        raise Exception('dir must be specified with analysis flag if fit flag is disabled')

    return opts


class SimulSingleCell( object):
    
    def __init__(self, opts):

        self.opts = opts 
        self.path_params = []
        self.path_analysis = []

        # filenames
        self.fname = {
                'initparams' : "initParams",
        }

        # sbatch options
        self.slurm= {
                'routine' : 'multicore', # singlecore or multicore
                'account' : 'ucb-summit-smr',
                'time' : '16:00:00',
                'qos' : 'condo',
                'partition' : 'shas',
                'jobname' : 'SingleCell',
                'output' : 'sim.log',
                'error' : 'sim.err',
                # Define architecture of clusters
                'coresPerNode' : 24,
                'socksPerNode' : 2,
                # 'username' : 'saan8193@colorado.edu',
                # 'mailtype' : 'FAIL',
        }

        # paths for different servers
        paths_workdir = {
                'Summit' : '/projects/saan8193/ImageAnalysis/SingleCell',
                'Local' : '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell'
        }
        paths_launch = {
                'Summit' : '/projects/saan8193',
                'Local' : '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell'
        }

        # set the correct working directory
        self.opts.loc = 'Summit'
        if self.opts.loc == 'Summit':
            self.workdir = paths_workdir['Summit']
            self.launchdir = paths_launch['Summit']
        elif self.opts.loc == 'Local':
            self.workdir = paths_workdir['Local']
            self.launchdir = paths_launch['Summit']
        elif self.opts.loc == 'Rumor':
            raise ValueError('RunSimulSingleCell: rumor not set up yet.')
        else:
            raise ValueError('RunSimulSingleCell: acceptable configuration input argument options are Summit, Local and Rumor.')

        os.chdir( self.workdir) 


    def Launch(self):
        # Launcher for multiple single cell fits

        # Prepare fits
        if self.opts.fit:
            # Initialize params for the different cells
            self.InitializeFit()
            
        # Prepare fits
        if self.opts.analyze:
            # Initialize params for the different cells
            self.InitializeAnalysis()

        # create a bash script for executing the jobs
        jobStrings = self.WriteLaunchScript()

        # Launch individually
        if (self.opts.fit or self.opts.analyze):
            for spath, jobString in zip( self.path_params, jobStrings):
                os.chdir( os.path.split(spath)[0])
                with open('jobscript.sh', 'w') as f:
                    f.write( jobString)
                if not self.opts.dryrun:
                    subprocess.call(["sbatch", "jobscript.sh"])
        # if not self.opts.dryrun and (self.opts.fit or self.opts.analyze):
           # status = call(['sbatch', self.fname['launch']])  


    def InitializeFit(self):
        # initialize fit : initialize params.mat for the different cells and get their path locations

        # check if initParams exists in the current directory
        if os.path.exists( self.fname['initparams']+'.m' ) == False:
            raise ImportError('RunSimulSingleCell: {0}.m does not exist in the current working directory'.format( self.fname['initparams'] ) )
        
        print( 'Initializing Fit: running initParams.m to save parameters for the segmented cells :')
        
        # run parameter initialization and get paths to saved parameter files
        if self.opts.loc == 'Summit':
            opts_params = { 'LOC': 'Summit', 'CFG': 'RELEASE'}
        elif self.opts.loc == 'Rumor':
            opts_params = { 'LOC': 'Rumor', 'CFG': 'RELEASE'}

        eng = matlab.engine.start_matlab()
        self.path_params = getattr( eng, self.fname['initparams'])( opts_params)
        for path in self.path_params:
            print('{0}'.format( path) ) 

        # number of jobs based on length of pathParams
        self.n_jobs = len( self.path_params)


    def InitializeAnalysis(self):
        # initialize analysis : get the paths of the cells to analyze

        print( 'Initializing Analysis : getting paths for analysis :')

        # If also launching fits, then just use the path to saved params to get the parent directory where analysis will be run
        if self.opts.fit:
            [self.path_analysis.append( os.path.dirname( ppath) ) for ppath in self.path_params]

        # Otherwise, if directory not specified, prompt an error.
        elif not self.opts.dir:
            raise ValueError('For analysis only without launching fits, user must specify directory string to match to')

        # Otherwise, use input argument to match folder names
        else:

            # Parent path where cell of interest is located 
            apath = os.path.dirname( self.opts.dir)
            apath = os.path.abspath( apath)

            # Get all folders in the parent directory of interested file
            alldirs = [os.path.join(apath, o) for o in os.listdir(apath) if os.path.isdir(os.path.join(apath,o))]

            # match the regex( possibly) analyzeDir name to folders that exist in apath
            self.path_analysis = [ cell for cell in alldirs if re.findall( self.opts.dir, cell) ];

            if not self.path_analysis:
                raise ValueError('There were no matching directories for analysis in directory specified in argument')

            # number of jobs based on length of analysis paths 
            self.n_jobs = len( self.path_analysis)

    def WriteLaunchScript(self):
        # create sbatch bash script for execution on summit
        
        jobStrings = []
        nTasks = 1
           
        # Define jobString
        jobStringDef = """#!/bin/bash

#SBATCH --job-name={0}
#SBATCH --qos={1}
#SBATCH --partition={2}
#SBATCH --account={3}
#SBATCH --output=fit.log
#SBATCH --error=fit.err
#SBATCH --time={4}
#SBATCH --nodes={5}
#SBATCH --ntasks={6}

export SCRATCH=/scratch/summit/saan8193
mkdir -p $SCRATCH/$SLURM_JOB_ID

module purge
module load matlab/R2018b

"""
        # Loop over seeds and make job strings to launch
        for spath in self.path_params:

            # Find number of nodes and number of processors/task
            if self.slurm['routine'] == 'singlecore':
                nCpuPerTask = 1
                nNodes = int( np.ceil(nTasks/(float(self.slurm['coresPerNode'])/nCpuPerTask)) )
            elif self.slurm['routine'] == 'multicore':
                nCpuPerTask = 24 
                nNodes = int( np.ceil(nTasks/(float(self.slurm['coresPerNode'])/nCpuPerTask)) )

            # Jobname : SimName_SeedNumber
            jobName = '__'.join( spath.split('/')[-2:] )

            # Write jobString 
            jobString = jobStringDef.format( self.slurm['jobname'], self.slurm['qos'], self.slurm['partition'], self.slurm['account'], self.slurm['time'], nNodes, nCpuPerTask) 

            jobString = jobString + 'matlab -nodesktop -r "clear all; setenv({3},{4}); cd {0}; addpath({1}); singleCell({2})"\n'.format( repr('/projects/saan8193/ImageAnalysis/SingleCell'), repr('classes'), repr(spath), repr('TZ'), repr('America/Denver'))
            jobString = jobString + '\nrm -rf $SCRATCH/$SLURM_JOB_ID\n'
            jobStrings += [jobString]

        return jobStrings

if __name__ == '__main__':
    
    # parse arguments
    opts = parse_args()
    # initialize multiple single cell jobs
    simulCells = SimulSingleCell(opts)
    # launch single cell jobs 
    simulCells.Launch()

