#!/usr/bin/env python

import os
import sys
import pdb
import matlab.engine
import argparse
import re
from subprocess import call

'''
Name : RunSimulSingleCell.py
Description : Runs multiple fits on Summit by distirbuting them across multiple cores and nodes as needed.
Usage : RunSimulSingleCell.py summit
'''

def parse_args():
    # parse comand line arguments

    parser = argparse.ArgumentParser( prog='SimulSingleCell.py', description='Launcher for multiple matlab single cell fits on summit') 
    
    # configuration of launcher and the matlab fit program
    parser.add_argument('--loc', type=str, default='Summit', 
            help='define the configuration environment for both the launcher and the matlab fit program')

    # Launch option: this must be enabled for launching runs
    parser.add_argument('-L', '--launch', action='store_true',  
            help='launch the jobs via sbatch')

    # Analyze flag: this must be enabled to run analysis on each cell after completion of fitting.
    parser.add_argument('-A', '--analyzeDir', type=str, nargs='?', default = 'TBD',
            help='Analyze the fitted cell data present in cells in specified directory.')

    opts = parser.parse_args()
    return opts


class SimulSingleCell( object):
    
    def __init__(self, opts):

        self.opts = opts 

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

        # initialize params for the different cells
        self.InitializeParams()

        # Find the number of cores and nodes these jobs should utilize
        self.FindNumberCoresNodes()

        # write list of bash commands to be executed by loadbalancer
        self.WriteJobsFile()

        # create a bash script for executing the jobs
        self.WriteLaunchScript()

        # Launch and Analyze the fits
        if self.opts.launch or self.opts.analyze:
           status = call(['sbatch', self.fname['launch']])  


    def InitializeParams(self):
        # initialize params for the different cells  and save their path locations

        # check if initParams exists in the current directory
        if os.path.exists( self.fname['initparams']+'.m' ) == False:
            raise ImportError('RunSimulSingleCell: {0}.m does not exist in the current working directory'.format( self.fname['initparams'] ) )
        
        print( 'Initalizing parameters for the segmented cells :')
        
        # run parameter initialization and get paths to saved parameter files
        opts_params = { 'LOC': 'Summit', 'CFG': 'RELEASE'}
        eng = matlab.engine.start_matlab()
        self.path_params = getattr( eng, self.fname['initparams'])( opts_params)
        for path in self.path_params:
            print('{0}'.format( path) ) 

        # number of jobs based on length of pathParams
        self.n_jobs = len( self.path_params)


    def FindNumberCoresNodes(self):
        # figure out number of cores and nodes based on simple division and modulus operation

        # figure out nodes and cores numbers
        self.n_nodes = 1 + self.n_jobs // self.n_cores_per_node
        self.n_cores = self.n_jobs % self.n_cores_per_node 

        print('Number of Jobs = {0}\nNumber of Nodes = {1}\nNumber of Cores = {2}'.format(self.n_jobs, self.n_nodes, self.n_cores) )

    def WriteJobsFile(self):
        # writes a list of bash commands for each Cell to be executed by loadbalancer

        # open file
        f_jobs = open(self.fname['jobs'], "w+")
        print( 'Writing jobs file : {0}'.format(self.fname['jobs']) )
        self.WriteLaunchAnalyze
        print( 'Successful write to jobs file!') 


    def WriteLaunchAnalyze(self):
        # write launch and subsequent analysis jobs information 

        # Loop over cells and echo command into jobs file
        for idx, ppath in enumerate( self.path_params):

            # Prefix for job command
            cmd_pre = 'matlab -nodesktop -r "clear all; addpath( genpath( "classes"));'

            # Suffix for job command
            cmd_post = '"\n'

            # Launch command
            cmd_launch = 'singleCell({0});'.format( repr(ppath) )

            # Analyze command 
            if self.opts.analyzeDir:
                self.GetPathsAnalysis
                cmd_analysis = 'AnalysisSingleCell.AnalyzeSingle({0});'.format( repr( self.path_analysis[idx] ) )

            # Write launch or analysis command
            if self.opts.launch and self.opts.analyzeDir:
                job_cmd = cmd_pre + cmd_launch + cmd_analysis + cmd_post

            elif self.opts.launch:
                job_cmd = cmd_pre + cmd_launch + cmd_post

            elif self.opts.analyze:
                job_cmd = cmd_pre + cmd_analysis + cmd_post

            # write this out and close the file
            f_jobs.write( job_cmd) 
        
        f_jobs.close()

    def GetPathsAnalysis(self):
        # use the analysisDir argument to get paths of cells to analyze if possible

        # If also launching fits, then just use the path to saved params to get the parent directory where analysis will be run
        if self.opts.launch:
            for idx, ppath in enumerate( self.path_params):
                self.path_analysis[idx] = os.path.dirname( ppath)


        # Otherwise, if directory not specified, prompt an error.
        elif not self.opts.analyzeDir:
            ValueError('For analysis only without launching fits, user must specify directory string to match to')

        # Otherwise, use input argument to match folder names
        else:

            # Parent path where cell of interest is located 
            apath = os.path.dirname( self.opts.analyzeDir)
            apath = os.path.abspath( apath)

            # Get all folders in the parent directory of interested file
            alldirs = [os.path.join(apath, o) for o in os.listdir(apath) if os.path.isdir(os.path.join(apath,o))]

            # match the regex( possibly) analyzeDir name to folders that exist in apath
            self.path_analysis = [ cell for cell in alldirs if re.findall( self.opts.analyzeDir, cell) ];

            if not self.path_analysis:
                ValueError('There were no matching directories for analysis in directory specified in argument')
                     
        print( self.path_analysis)
                    

    def WriteLaunchScript(self):
        # create sbatch bash script for execution on summit
        
        # open file
        f_launch = open(self.fname['launch'], "w+")
        print( 'Writing launch file : {0}'.format(self.fname['launch']) )

        # write sbatch options
        self.WriteSbatchHeader

        # job execution command
        f_launch.write( '\n# launch cell fitting jobs\n')
        f_launch.write( 'mpirun lb {0}'.format( self.fname['jobs']) )

        # close the file
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

        # modules to load
        modules = ['intel', 'impi', 'python', 'loadbalance', 'matlab']

        # load modules
        f_launch.write( '\n# clearing all modules\n')
        f_launch.write( 'module purge\n')
        f_launch.write( '\n# load required modules\n')
        for module in modules:
            f_launch.write( 'module load {0}\n'.format(module) )

        # close the file
        f_launch.close()


if __name__ == '__main__':
    
    # parse arguments
    opts = parse_args()

    # initialize multiple single cell jobs
    simulCells = SimulSingleCell(opts)

    # launch single cell jobs 
    simulCells.Launch()

