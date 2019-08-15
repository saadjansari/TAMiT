SCIANT : Single Cell Image Analysis Tool

Summit:

1. Movies (.nd2) need to be pre-segmented locally. Use segmentMovie() function. Segmented cells will be saved as .mat files
2. Segmented cells need to be transfered to Summit (Remote Server)
    scp -r <pathToSegmentedCell> <user>@login.rc.colorado.edu:/projects/<user>/<pathToCells>
    example:
        scp -r 998_150msR_100G_trig_7Z_001_cells saan8193@login.rc.colorado.edu:"/projects/saan8193/ImageAnalysis/FY\\ Datasets/."
3. Log in to summit (remote server) and ssh into a compile node
    ssh <user>@login.rc.colorado.edu
    ssh shas0137
4. Navigate to 



