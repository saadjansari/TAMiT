function simAndFitBatch(spath)

    addpath('classes');
    addpath(genpath('functions'))
    if nargin == 0
        spath = [];
    end
    
    nTrials = 200;
    snrVals = [2,4,8,16,32];
    test_type = 'Monopolar';
    savepath = [spath, filesep, 'simFitdata_', test_type, '.mat'];

    launch_summit = 1;
    pc = parcluster('local');
    if launch_summit
        % explicitly set the JobStorageLocation to a temp directory
        pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_JOB_ID'));
        parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))
    else
        parpool(pc);
    end

    dxyz = cell( length(snrVals), nTrials);
    spb_xyz = cell( length(snrVals), nTrials);
    FP = cell( length(snrVals), nTrials);
    TN = cell( length(snrVals), nTrials);
    lenTN = cell( length(snrVals), nTrials);
    lenTP = cell( length(snrVals), nTrials);
    lenFP = cell( length(snrVals), nTrials);

    for nidx = 1 : length( snrVals)
        parfor jT = 1 : nTrials

            dat = simAndFit( snrVals(nidx), test_type);
            dxyz{nidx, jT} = dat.dxyz;
            spb_xyz{nidx, jT} = dat.spb_xyz;
            FP{nidx,jT} = dat.FP;
            TN{nidx,jT} = dat.TN;
            lenTN{nidx, jT} = dat.lenTN;
            lenTP{nidx, jT} = dat.lenTP;
            lenFP{nidx, jT} = dat.lenFP;


            fid = fopen([spath, filesep, 'running.txt'],'w');
            fprintf(fid,['noise = ' num2str(snrVals(nidx)) '\nTrial = ' num2str(jT) '\n']);
            fclose(fid);

        end
    end

    % combine data and save
    save(savepath, 'nTrials', 'snrVals', 'dxyz', 'spb_xyz', 'FP', 'TN', 'lenTN', 'lenTP', 'lenFP')

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