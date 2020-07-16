function dat = simAndFit(noise,type)

    addpath('classes');
    addpath(genpath('functions'))

    %% simulate an image
    [imNoisy, mainObj, params] = simCell(noise, type);
    
    
    %% Fit
    disp('Estimating features...')
    % Estimate the features based on an estimation routine (defined in specialized sub-class )
    switch params.type
        case 'Monopolar'
            
            costAccept = 2; 
            
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
            fpar.fitExploreSpeed = 1;
            fpar.fit2DFirst = 0;
            fpar.channel = 1; fpar.time = 1; fpar.saveDirectory = [pwd, filesep, 'simdata']; fpar.timeReversal = 0;
            feature.updateSubFeatures(); feature.syncFeaturesWithMap(); 
            
            % find best length
            for jf = 2 : length(feature.featureList{1}.featureList)
                L0 = feature.featureList{1}.featureList{jf}.length;
                Ls = linspace( L0-10, L0+10, 20); Ls( Ls < 5) = [];
                resids = [];
                for el = Ls
                    feature.featureList{1}.featureList{jf}.SetLength(el);
                    ims = feature.simulateAll( imNoisy, feature.ID);
                    resids = [resids, sum( (ims(:) - imNoisy(:) ).^2)];
                end
                [~, idx] = min(resids); feature.featureList{1}.featureList{jf}.SetLength( Ls(idx) );
            end
            
            % delete microtubules less than length = 7
            rmLine = [];
            for jf = 2 : length(feature.featureList{1}.featureList)
                if feature.featureList{1}.featureList{jf}.length < 7
                    rmLine = [ rmLine jf];
                end
            end
            feature.featureList{1}.featureList(rmLine) = []; feature.featureList{1}.numFeatures = length( feature.featureList{1}.featureList);
            feature.syncFeaturesWithMap(); 
%             figure; imagesc( max(imNoisy,[],3) ); axis equal; colormap gray; hold on; feature.displayFeature( gca)

            % Fit Features via Fit Engine.
            fitEngine = FitEngine( imNoisy, feature, fpar);
            fitEngine = fitEngine.Optimize();
            finalFit = fitEngine.GetFeature();

            % Basic feature linking
            
            % SPB data
            dat.spb_xyz = norm( finalFit.featureList{1}.featureList{1}.position - mainObj.featureList{1}.featureList{1}.position);
            dat.spb_true = mainObj.featureList{1}.featureList{1}.position;
            
            % Line data
            dat.noise = noise; 
            costs = zeros( length(finalFit.featureList{1}.featureList)-1, length(mainObj.featureList{1}.featureList)-1);
            for j1 = 1 : size( costs, 1)
                for j2 = 1 : size( costs, 2)
                    e1 = finalFit.featureList{1}.featureList{1+j1}.endPosition;
                    e2 = mainObj.featureList{1}.featureList{1+j2}.endPosition;
                    costs(j1,j2) = norm( e1-e2);
                end
            end
            
            % for features found, how good are they
            dxyz = []; dat.FP = 0; jList = []; dat.lenTP = [];
            for j1 = 1 : size( costs, 1)
                [minn, midx] = min( costs(j1,:) );
                if minn < costAccept
                    jList = [jList midx];
                    diff = finalFit.featureList{1}.featureList{1+j1}.endPosition - mainObj.featureList{1}.featureList{1+midx}.endPosition;
                    dxyz = [ dxyz ; abs(diff)];
                    dat.lenTP = [ dat.lenTP, norm( [0.1 0.1 0.5].*(finalFit.featureList{1}.featureList{1+j1}.endPosition - finalFit.featureList{1}.featureList{1+j1}.startPosition))];
                else
                    dat.FP = dat.FP + 1;
                    dat.lenFP = [ dat.lenFP, norm( [0.1 0.1 0.5].*(finalFit.featureList{1}.featureList{1+j1}.endPosition - finalFit.featureList{1}.featureList{1+j1}.startPosition))];
                end
            end
            dat.TN = size(costs,2) - size(dxyz,1);
            dat.dxyz = dxyz;
            
            % Line data
            dat.negAmp = []; dat.lenTN = [];
            for j2 = 1 : size( costs, 2)
                if ~any(jList == j2)
                    dat.negAmp = [dat.negAmp, mainObj.featureList{1}.featureList{1+j2}.amplitude];
                    dat.lenTN = [ dat.lenTN, norm( [0.1 0.1 0.5].*(mainObj.featureList{1}.featureList{1+j1}.endPosition - mainObj.featureList{1}.featureList{1+j1}.startPosition))];
                end
            end
            
        case 'Spindle'
            
            
        case 'IMTBank'
            
    end
%     figure;
%     subplot(121); imagesc( max(imNoisy,[],3)); axis equal; colormap gray; title(['Noise = ' num2str(noise)]); xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     subplot(122); imagesc( max( finalFit.simulateAll(imNoisy, finalFit.ID),[],3)); axis equal; colormap gray; title('Simulated'); xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     drawnow; pause(1)
    
end