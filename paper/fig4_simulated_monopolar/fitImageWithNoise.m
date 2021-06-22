function finalFit = fitImageWithNoise(imSim, mainObj, params, snr)

    addpath(genpath('../../classes'));
    addpath(genpath('../../functions'))

    %% add noise to image
    imNoisy = mat2gray( add_poisson_noise( imSim, snr) );
    
    %% Fit
    disp('Estimating features...')
    % Estimate the features based on an estimation routine (defined in specialized sub-class )

    costAccept = 4;

    par.display = 1;
    par.pole =1;
    par.mt =1;
    props = Cell.GetFeatureProps();
    feature = MonopolarCell.findFeaturesDeNovo_MT( imNoisy, par, props.monopolarAster);

    fpar.runFit = 1;
    fpar.runLocalOptimization = 0;
    fpar.runGlobalOptimization = 1;
    fpar.runFeatureNumberOptimization = 1;
    fpar.useParallel = true;
    fpar.state = 'RELEASE';
    fpar.display = 0;
    fpar.alpha = 0.4;
    fpar.fitExploreSpeed = 0;
    fpar.fit2DFirst = 1;
    fpar.channel = 1; fpar.time = 1; fpar.saveDirectory = [pwd, filesep, 'simdata']; fpar.timeReversal = 0;
    fpar.scaleParameters = 1;
    global COUNTER
    COUNTER = 1;
    feature.ID = COUNTER;
    COUNTER = COUNTER + 1;
    feature.syncFeatures();
    feature.preOptimize();

    % Fit Features via Fit Engine.
    fitEngine = FitEngine( imNoisy, feature, fpar);
    fitEngine = fitEngine.Optimize();
    finalFit = fitEngine.GetFeature();
    
    % Basic feature linking

%     % SPB data
%     dat.spb_xyz = norm( finalFit.featureList{1}.featureList{1}.position - mainObj.featureList{1}.featureList{1}.position);
%     dat.spb_true = mainObj.featureList{1}.featureList{1}.position;
% 
%     % Line data
%     dat.snr = snr; 
%     costs = zeros( length(finalFit.featureList{1}.featureList)-1, length(mainObj.featureList{1}.featureList)-1);
%     for j1 = 1 : size( costs, 1)
%         for j2 = 1 : size( costs, 2)
%             e1 = finalFit.featureList{1}.featureList{1+j1}.endPosition;
%             e2 = mainObj.featureList{1}.featureList{1+j2}.endPosition;
%             costs(j1,j2) = norm( e1-e2);
%         end
%     end
% 
%     % for features found, how good are they
%     dxyz = []; dat.FP = 0; jList = []; dat.lenTP = []; dat.lenFP = [];
%     for j1 = 1 : size( costs, 1)
%         [minn, midx] = min( costs(j1,:) );
%         if minn < costAccept
%             jList = [jList midx];
%             diff = finalFit.featureList{1}.featureList{1+j1}.endPosition - mainObj.featureList{1}.featureList{1+midx}.endPosition;
%             dxyz = [ dxyz ; abs(diff)];
%             dat.lenTP = [ dat.lenTP; finalFit.featureList{1}.featureList{1+j1}.endPosition - finalFit.featureList{1}.featureList{1+j1}.startPosition];
%         else
%             dat.FP = dat.FP + 1;
%             dat.lenFP = [ dat.lenFP; finalFit.featureList{1}.featureList{1+j1}.endPosition - finalFit.featureList{1}.featureList{1+j1}.startPosition];
%         end
%     end
%     dat.TN = size(costs,2) - size(dxyz,1);
%     dat.dxyz = dxyz;
% 
%     % Line data
%     dat.lenTN = [];
%     for j2 = 1 : size( costs, 2)
%         if ~any(jList == j2)
%             dat.lenTN = [ dat.lenTN; mainObj.featureList{1}.featureList{1+j2}.endPosition - mainObj.featureList{1}.featureList{1+j2}.startPosition];
%         end
%     end
        
%     figure;
%     subplot(121); imagesc( max(imNoisy,[],3)); axis equal; colormap gray; title(['Noise = ' num2str(noise)]); xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     subplot(122); imagesc( max( finalFit.simulateAll(imNoisy, finalFit.ID),[],3)); axis equal; colormap gray; title('Simulated'); xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     drawnow; pause(1)
    
end
