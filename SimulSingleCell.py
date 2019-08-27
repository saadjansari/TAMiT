#!/usr/bin/env python

import os
import sys
import pdb
import matlab.engine
import argparse
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
    parser.add_argument('-cfg', '--config', type=str, default='Summit', 
            help='define the configuration environment for both the launcher and the matlab fit program')

    # configuration of launcher and the matlab fit program
    parser.add_argument('-L', '--launch', action='store_true',  
            help='launch the jobs via sbatch')

    opts = parser.parse_args()
    return opts


class SimulSingleCell( object):
    
    def __init__(self, opts):

        self.opts = opts 

        # number of jobs, nodes and cores for this job
        self.n_nodes = []
        self.n_cores = []
        self.n_jobs = []
        self.n_cores_per_node = 12

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
                'time' : '05:00:00',
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
        if self.opts.config == 'Summit':
            self.workdir = paths_workdir['Summit']
            self.launchdir = paths_launch['Summit']
            
        elif self.opts.config == 'Local':
            self.workdir = paths_workdir['Local']
            self.launchdir = paths_launch['Summit']
            
        elif self.opts.config == 'Rumor':
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

        # Launch the fits
        if self.opts.launch:
            print( 'Launching simulataneous cell fits...')
            # os.chdir( self.launchdir) 
            status = call(['sbatch', self.fname['launch']])  

    def InitializeParams(self):
        # initialize params for the different cells  and save their path locations

        # check if initParams exists in the current directory
        if os.path.exists( self.fname['initparams']+'.m' ) == False:
            raise ImportError('RunSimulSingleCell: {0}.m does not exist in the current working directory'.format( self.fname['initparams'] ) )
        
        print( 'Initalizing parameters for the segmented cells :')
        
        # run parameter initialization and get paths to saved parameter files
        eng = matlab.engine.start_matlab()
        self.path_params = getattr( eng, self.fname['initparams'])( self.opts.config )
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


    def WriteLaunchScript(self):
        # create sbatch bash script for execution on summit
        
        # open file
        f_launch = open(self.fname['launch'], "w+")
        print( 'Writing launch file : {0}'.format(self.fname['launch']) )

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

        f_launch.write( '\n# clearing all modules\n')
        f_launch.write( 'module purge\n')
        f_launch.write( '\n# load required modules\n')
        for module in modules:
            f_launch.write( 'module load {0}\n'.format(module) )

        # job execution command
        f_launch.write( '\n# launch cell fitting jobs\n')
        f_launch.write( 'mpirun lb {0}'.format( self.fname['jobs']) )

        # close the file
        f_launch.close()
        print( 'Successful write to launch file!') 


    def WriteJobsFile(self):
        # writes a list of bash commands for each Cell to be executed by loadbalancer

        # set working directory
        # os.chdir( '/projects/saan8193/')

        # open file
        f_jobs = open(self.fname['jobs'], "w+")
        print( 'Writing jobs file : {0}'.format(self.fname['jobs']) )

        # Loop over cells and echo command into jobs file
        for path in self.path_params:

            # command for launching a single job
            job_cmd = 'matlab -nodesktop -r "clear all; singleCell({0});"\n'.format( repr(path) )

            # write this out and close the file
            f_jobs.write( job_cmd) 
        
        f_jobs.close()
        print( 'Successful write to jobs file!') 


if __name__ == '__main__':
    
    # parse arguments
    opts = parse_args()

    # initialize multiple single cell jobs
    simulCells = SimulSingleCell(opts)

    # launch single cell jobs 
    simulCells.Launch()



