#!/usr/bin/env python

import os
import sys
import pdb
import matlab.engine
import argparse
import re
import pdb
from subprocess import call

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

        # number of jobs, nodes and cores for this job
        self.n_nodes = []
        self.n_cores = []
        self.n_jobs = []
        self.n_cores_per_node = 24

        # filenames
        self.fname = {
                'initparams' : "initParams",
                'jobs' : "simul_single_cell_jobs",
                'launch' : "simul_single_cell_launch.sh"
        }

        # sbatch options
        self.sbatch = {
                'username' : 'saan8193@colorado.edu',
                'mailtype' : 'FAIL',
                'account' : 'ucb-summit-smr',
                'jobname' : 'SimulSingleCell',
                'qos' : 'condo',
                'partition' : 'shas',
                'time' : '10:00:00',
                'output' : 'simul_single_cell.out'
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

        # Find the number of cores and nodes these jobs should utilize
        self.FindNumberCoresNodes()

        # write list of bash commands to be executed by loadbalancer
        self.WriteJobsFile()

        # create a bash script for executing the jobs
        self.WriteLaunchScript()

        # Launch and Analyze the fits
        if not self.opts.dryrun and (self.opts.fit or self.opts.analyze):
           status = call(['sbatch', self.fname['launch']])  


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


    def FindNumberCoresNodes(self):
        # figure out number of cores and nodes based on simple division and modulus operation

        # figure out nodes and cores numbers
        self.n_nodes = 1 + self.n_jobs // self.n_cores_per_node
        self.n_cores = self.n_jobs % self.n_cores_per_node 

        print('Number of Jobs = {0}\nNumber of Nodes = {1}\nNumber of Cores = {2}'.format(self.n_jobs, self.n_nodes, self.n_cores) )

    def WriteJobsFile(self):
        # writes a list of bash commands for each Cell to be executed by loadbalancer

        print( 'Writing jobs file : {0}'.format(self.fname['jobs']) )

        # open file
        f_jobs = open(self.fname['jobs'], "w+")

        # Prefix for job command
        cmd_pre = 'matlab -nodesktop -r "clear all; addpath( genpath( "classes"));'

        # Suffix for job command
        cmd_post = '"\n'

        # Fit and Analyze
        if self.opts.fit and self.opts.analyze:
            
            for pathp, patha in zip(self.path_params, self.path_analysis):

                cmd_launch = 'singleCell({0});'.format( repr(pathp) )
                cmd_analysis = 'AnalysisSingleCell.AnalyzeSingle({0});'.format( repr( patha ) )
                job_cmd = cmd_pre + cmd_launch + cmd_analysis + cmd_post
                f_jobs.write( job_cmd) 
        
        # Fit Only
        elif self.opts.fit:
            
            for pathp in self.path_params:

                cmd_launch = 'singleCell({0});'.format( repr(pathp) )
                job_cmd = cmd_pre + cmd_launch + cmd_post
                f_jobs.write( job_cmd) 
        
        # Analyze Only
        elif self.opts.analyze:

            for patha in self.path_analysis:

                cmd_analysis = 'AnalysisSingleCell.AnalyzeSingle({0});'.format( repr( patha ) )
                job_cmd = cmd_pre + cmd_analysis + cmd_post
                f_jobs.write( job_cmd) 

        f_jobs.close()
        print( 'Successful write to jobs file!') 


    def WriteLaunchScript(self):
        # create sbatch bash script for execution on summit
        
        print( 'Writing launch file : {0}'.format(self.fname['launch']) )

        # Write sbatch options
        self.WriteSbatchHeader()

        f_launch = open(self.fname['launch'], "a")
        # Load modules
        modules = ['intel', 'impi', 'python', 'loadbalance', 'matlab']
        f_launch.write( '\n# clearing all modules\n')
        f_launch.write( 'module purge\n')
        f_launch.write( '\n# load required modules\n')
        for module in modules:
            f_launch.write( 'module load {0}\n'.format(module) )

        # Job execution command
        f_launch.write( '\n# launch cell fitting jobs\n')
        f_launch.write( 'mpirun lb {0}'.format( self.fname['jobs']) )

        # Close the file
        f_launch.close()
        print( 'Successful write to launch file!') 


    def WriteSbatchHeader(self):
        # create sbatch bash script for execution on summit
        
        # open file
        f_launch = open(self.fname['launch'], "w+")

        # identiy bash evironment and define sbatch options
        f_launch.write( '#!/bin/bash\n\n')
        f_launch.write( '#SBATCH --job-name {0}\n'.format( self.sbatch['jobname']) )
        f_launch.write( '#SBATCH --qos={0}\n'.format( self.sbatch['qos']) )
        f_launch.write( '#SBATCH --partition={0}\n'.format( self.sbatch['partition']) )
        f_launch.write( '#SBATCH --mail-type={0}\n'.format( self.sbatch['mailtype']) )
        f_launch.write( '#SBATCH --mail-user={0}\n'.format( self.sbatch['username']) )
        f_launch.write( '#SBATCH --account={0}\n'.format( self.sbatch['account']) )
        f_launch.write( '#SBATCH --time {0}\n'.format( self.sbatch['time']) )
        f_launch.write( '#SBATCH --output {0}\n'.format( self.sbatch['output']) )
        f_launch.write( '#SBATCH --nodes {0}\n'.format(self.n_nodes) )
        f_launch.write( '#SBATCH --ntasks {0}\n'.format(self.n_jobs) )

        # close the file
        f_launch.close()


if __name__ == '__main__':
    
    # parse arguments
    opts = parse_args()

    # initialize multiple single cell jobs
    simulCells = SimulSingleCell(opts)

    # launch single cell jobs 
    simulCells.Launch()

