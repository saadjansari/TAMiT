function dat = simAndFit(snr,type)

    addpath('classes');
    addpath(genpath('functions'))

    %% simulate an image
    [imNoisy, mainObj, params] = simCell(snr, type);
    imNoisy = mat2gray(imNoisy);
    %% Fit
    disp('Estimating features...')
    % Estimate the features based on an estimation routine (defined in specialized sub-class )
    switch params.type
        case 'Monopolar'
            
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
            fpar.fitExploreSpeed = 1;
            fpar.fit2DFirst = 1;
            fpar.channel = 1; fpar.time = 1; fpar.saveDirectory = [pwd, filesep, 'simdata']; fpar.timeReversal = 0;
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
            
            % SPB data
            dat.spb_xyz = norm( finalFit.featureList{1}.featureList{1}.position - mainObj.featureList{1}.featureList{1}.position);
            dat.spb_true = mainObj.featureList{1}.featureList{1}.position;
            
            % Line data
            dat.snr = snr; 
            costs = zeros( length(finalFit.featureList{1}.featureList)-1, length(mainObj.featureList{1}.featureList)-1);
            for j1 = 1 : size( costs, 1)
                for j2 = 1 : size( costs, 2)
                    e1 = finalFit.featureList{1}.featureList{1+j1}.endPosition;
                    e2 = mainObj.featureList{1}.featureList{1+j2}.endPosition;
                    costs(j1,j2) = norm( e1-e2);
                end
            end
            
            % for features found, how good are they
            dxyz = []; dat.FP = 0; jList = []; dat.lenTP = []; dat.lenFP = [];
            for j1 = 1 : size( costs, 1)
                [minn, midx] = min( costs(j1,:) );
                if minn < costAccept
                    jList = [jList midx];
                    diff = finalFit.featureList{1}.featureList{1+j1}.endPosition - mainObj.featureList{1}.featureList{1+midx}.endPosition;
                    dxyz = [ dxyz ; abs(diff)];
                    dat.lenTP = [ dat.lenTP; finalFit.featureList{1}.featureList{1+j1}.endPosition - finalFit.featureList{1}.featureList{1+j1}.startPosition];
                else
                    dat.FP = dat.FP + 1;
                    dat.lenFP = [ dat.lenFP; finalFit.featureList{1}.featureList{1+j1}.endPosition - finalFit.featureList{1}.featureList{1+j1}.startPosition];
                end
            end
            dat.TN = size(costs,2) - size(dxyz,1);
            dat.dxyz = dxyz;
            
            % Line data
            dat.lenTN = [];
            for j2 = 1 : size( costs, 2)
                if ~any(jList == j2)
                    dat.lenTN = [ dat.lenTN; mainObj.featureList{1}.featureList{1+j2}.endPosition - mainObj.featureList{1}.featureList{1+j2}.startPosition];
                end
            end
            
        case 'MitosisBud'
                        
            costAccept = 5;
            
            par.display = 1;
            par.spindleMT =1; par.astralMT=1; par.spindlePoles=1;
            par.mt =1;
            props = Cell.GetFeatureProps();
            feature = MitoticCellBud.findFeaturesDeNovo_MT( imNoisy, par, props.spindle);
            
            fpar.runFit = 1;
            fpar.runLocalOptimization = 0;
            fpar.runGlobalOptimization = 1;
            fpar.runFeatureNumberOptimization = 1;
            fpar.useParallel = true;
            fpar.state = 'RELEASE';
            fpar.display = 0;
            fpar.alpha = 0.1;
            fpar.fitExploreSpeed = 1;
            fpar.fit2DFirst = 1;
            fpar.channel = 1; fpar.time = 1; fpar.saveDirectory = [pwd, filesep, 'simdata']; fpar.timeReversal = 0;
            global COUNTER
            COUNTER = 1;
            feature.ID = COUNTER; 
            COUNTER = COUNTER + 1;
            feature.syncFeatures();
            feature.preOptimize();
            
%             figure; imagesc( max(imNoisy,[],3) ); axis equal; colormap gray; hold on; feature.displayFeature( gca)

            % Fit Features via Fit Engine.
            fitEngine = FitEngine( imNoisy, feature, fpar);
            fitEngine = fitEngine.Optimize();
            finalFit = fitEngine.GetFeature();
            
            % First match spbs
            costs = zeros(2);
            for j1 = 1 : size( costs, 1)
                for j2 = 1 : size( costs, 2)
                    e1 = finalFit.featureList{j1+1}.featureList{1}.position;
                    e2 = mainObj.featureList{j2+1}.featureList{1}.position;
                    costs(j1,j2) = norm( e1-e2);
                end
            end
            whichspb = zeros(1,2);
            [~,whichspb(1)] = min( costs(1, :)); [~,whichspb(2)] = min( costs(2, :));
            
            dat.snr = snr;
            % link each curved microtubule in each aster
            dat.dxyz = []; dat.FP = 0; dat.TN=0; jList = []; dat.lenTP = []; dat.lenFP = []; dat.lenTN = [];
            for jA = 1:2
                a1 = finalFit.featureList{1+jA};
                a2 = mainObj.featureList{1+whichspb(jA)};
                costs = zeros( length(a1.featureList)-1, length(a2.featureList)-1);
                for j1 = 1 : size( costs, 1)
                    for j2 = 1 : size( costs, 2)
                        c1 = a1.featureList{1+j1}.GetCoords(); c2 = a2.featureList{1+j2}.GetCoords();
                        e1 = c1(:,end)';
                        e2 = c2(:,end)';
                        costs(j1,j2) = norm( e1-e2);
                    end
                end
                dxyzt = [];
                % for features found, how good are they
                for j1 = 1 : size( costs, 1)
                    [minn, midx] = min( costs(j1,:) );
                    if minn < costAccept
                        jList = [jList midx];
                        c1 = a1.featureList{1+j1}.GetCoords(); c2 = a2.featureList{1+midx}.GetCoords();
                        diff = c1(:,end)'-c2(:,end)';
                        dxyzt = [ dxyzt ; abs(diff)];
                        dat.lenTP = [ dat.lenTP; a1.featureList{1+j1}.L];
                    else
                        dat.FP = dat.FP + 1;
                        dat.lenFP = [ dat.lenFP; a1.featureList{1+j1}.L];
                    end
                end
                dat.TN = dat.TN+size(costs,2) - size(dxyzt,1);
                dat.dxyz = [dat.dxyz; dxyzt];
            
                % Line data
                for j2 = 1 : size( costs, 2)
                    if ~any(jList == j2)
                        dat.lenTN = [ dat.lenTN; a2.featureList{1+j2}.L];
                    end
                end
            end
            
        case 'IMTBank'
            
    end
%     figure;
%     subplot(121); imagesc( max(imNoisy,[],3)); axis equal; colormap gray; title(['Noise = ' num2str(noise)]); xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     subplot(122); imagesc( max( finalFit.simulateAll(imNoisy, finalFit.ID),[],3)); axis equal; colormap gray; title('Simulated'); xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     drawnow; pause(1)
    
end
