runLocal = 0; % 1 for local, 0 for remote (i.e. Summit Supercomputer)

% define paths for running the code, loading the file, and saving all information
if runLocal == 0 
    runpath = '/projects/saan8193/Curved_Microtubule_Detector';
    filepath = '/projects/saan8193/FY Datasets/';
    savepath = '/scratch/summit/saan8193/Curved_Microtubule_Detector/results';
elseif runLocal == 1
    runpath = '/Users/saadjansari/Documents/Projects/Curved_Microtubule_Detector';
    filepath = '/Users/saadjansari/Documents/Projects/FY Datasets/';
    savepath = '/Users/saadjansari/Documents/Projects/Curved_Microtubule_Detector/results';
end


