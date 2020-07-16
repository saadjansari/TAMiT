function simAndFitBatch()

    addpath('classes');
    addpath(genpath('functions'))
    
    launch_summit = 1;
    pc = parcluster('local');
    if launch_summit
        % explicitly set the JobStorageLocation to a temp directory
        pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_JOB_ID'));
        parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))
    else
        parpool(pc);
    end
        
    nTrials = 1000;
    noiseVals = [0 0.2 0.4 0.6 0.8];
    
    dxyz = cell( length(noiseVals), nTrials);
    spb_xyz = cell( length(noiseVals), nTrials);
    FP = cell( length(noiseVals), nTrials);
    TN = cell( length(noiseVals), nTrials);
    lenTN = cell( length(noiseVals), nTrials);
    lenTP = cell( length(noiseVals), nTrials);
    lenFP = cell( length(noiseVals), nTrials);
    
    for nidx = 1 : length( noiseVals)
        parfor jT = 1 : nTrials
            
            dat = simAndFit( noiseVals(nidx), 'Monopolar');
            dxyz{nidx, jT} = dat.dxyz;
            spb_xyz{nidx, jT} = dat.spb_xyz;
            FP{nidx,jT} = dat.FP;
            TN{nidx,jT} = dat.TN;
            lenTN{nidx, jT} = dat.lenTN;
            lenTP{nidx, jT} = dat.lenTP;
            lenFP{nidx, jT} = dat.lenFP;
            
        end
    end
    
    % combine data and save
    save('simFitdata.mat', 'nTrials', 'noiseVals', 'dxyz', 'spb_xyz', 'FP', 'TN', 'lenTN', 'lenTP', 'lenFP')
    
%     figure;
%     dd = [];
%     for ii = 1:size( datt.dxyz,1)
%         dd = [dd, norm(datt.dxyz(ii,:))];
%     end
%     histogram(dd,15); xlabel('Error (voxel units)'); ylabel('Count'); title( 'Error of true positives (noise=0.4)')
%     
    
%     disp(['Sample size: ', num2str( length(dd) + datt.FP+datt.TN)])
%     disp(['True Positives: ', num2str( length(dd) )])
%     disp(['False Positives: ', num2str( datt.FP )])
%     disp(['True Negatives: ', num2str( datt.TN )])
    
end