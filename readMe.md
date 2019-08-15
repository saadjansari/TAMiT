SCIANT : Single Cell Image Analysis Tool

Summit:

1. Movies (.nd2) need to be pre-segmented locally. Use segmentMovie() function. Segmented cells will be saved as .mat files
2. Segmented cells need to be transfered to Summit (Remote Server)
    scp -r <pathToSegmentedCell> <user>@login.rc.colorado.edu:/projects/<user>/sciant/cells/.
3. Log in to summit (remote server)
    ssh <user>@login.rc.colorado.edu
4. Navigate to scratch sciant folder
    cd /scratch/summit/<user>/sciant/



