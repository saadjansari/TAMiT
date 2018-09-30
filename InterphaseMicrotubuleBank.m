classdef InterphaseMicrotubuleBank
% defines an organizing center for interphase microtubules. The bank can be thought of as a set of rules that governs the interactions between individual features. The bank is sort of an administrator at a given time, that figures out how individual features act when viewed together.
properties
    sourceImage
    mask
    dim
    sourceInfo
    numberOfMicrotubules
    display
    featureBank
    deadMicrotubules

end

% Main Methods
methods

    % Initialization
    % InterphaseMicrotubuleBank {{{
    function obj = InterphaseMicrotubuleBank( image3D, mask, cellNumber, currentTime, params)
    
        % All that is needed to initialize a bank is a source image( could be 2D or 3D).         
        dim = length( size(image3D));
        giveError( {dim, dim}, {'~=', '~='}, {2,3}, 'and', 'InterphaseMicrotubuleBank: 1st argument must be a 2D or a 3D image.')
        
        disp( 'Interphase Microtubule Bank created. ' )
        obj.sourceImage = image3D;
        obj.mask = mask;
        obj.dim = length( size(image3D) );
        obj.sourceInfo.fileName = params.meta.fileName;
        obj.sourceInfo.cellNumber = cellNumber;
        obj.sourceInfo.time = currentTime;
        obj.sourceInfo.cellCentroidInMovie = params.movieMatfile.cellCentroids( cellNumber, :);
        obj.sourceInfo.cellLocationImage = params.movieMatfile.imageNumbered;

    end
    % }}}

    % Estimation
    % EstimateMicrotubulesOld {{{
    function obj = EstimateMicrotubulesOld(obj, estParams)
        % EstimateMicrotubulePositions: estimates mt positions in a 2D/3D image.

        if obj.dim == 2               
            disp( 'Estimating microtubule locations...' )
        elseif obj.dim == 3
            disp( 'Estimating microtubule locations using a 2D MIP in the 3rd-dimension...' )
        end
        
        % extract the full cell image, and the segmented cell image
        imRegion = max( obj.sourceImage, [], 3);
        imCell = max( obj.sourceImage .* obj.mask, [], 3);
        imMask = logical( max( obj.mask, [], 3) );

        % Image Editing
        [imNetwork, imageGaussian] = EditMicrotubuleImageForEstimation( obj, imRegion, imMask);

        % Parse the image to separate out microtubules
        imMicrotubulesCell = ParseMicrotubulesFromNetwork( obj, imNetwork, imageGaussian, imMask, estParams);
         
        % find microtubule nucleation point by looking for the brightest intensity point along the microtubule.
        mtoc = findMTOC( obj, imMicrotubulesCell, imageGaussian, estParams);
        obj = obj;
        return
        % To estimate the lines, we'll take the following approach:
        % 1) Start at the MTOC (allow 2 nucleations on either side)
        % 2) For each nucleation, we'll find the angle that gives the maximum in a
        %    radially integrated signal, up to some max radius( which we will call
        %    the visibility)                                                                                                                                                        
        % 3) We'll propagate from the MTOC to a new point along the optimum angle,
        %    which a defined distance away from the MTOC (this we will call the
        %    step size). Throughout this process, we might look for optimum angles
        %    within a certain range related to the orientation of the Tube (this we
        %    will refer to as the field of vision).
        % 4) We will continue propagating until we can't find an optimum angle
        %    anymore (all our surroudings appear the same. This will be based on
        %    some optimality condition. At this point we will change our step size
        %    to half its original value and see if we can propogate a smaller
        %    distance forward.
        % 5) This will conclude our initial estimation phase.
        % 6) We will fit a polynomial through each initial estimate (quadratic or
        %    cubic). This will be our final estimate for the microtubules.
        
        mtBank = initializeMicrotubules( imMicrotubulesCell, imCell, mtoc)        

        % parameters
        polyOrder = 3;
        step = 2; vis = 8; fov = 60; % Step Size, Visibility, FieldOfVision
        
        % Estimate MT locations
        for jTube = 1 : numTubes
            mtBank(jTube) = Microtubule.EstimateMicrotubulePoints(mtBank( jTube),step,vis,fov);
            %     plotEstimatedPointsAll( mtBank( jTube) )
        end

        % Estimate polynomial curve coefficients
        for jTube = 1 : numTubes
            mtBank( jTube) = Microtubule.EstimateMicrotubuleCurve( mtBank( jTube), polyOrder);
            % PlotMicrotubuleCurve( mtBank( jTube), 'estimate')
        end
        
        % Estimate gaussian parameters
        for jTube = 1 : numTubes
            mtBank( jTube) = Microtubule.EstimateGaussianParameters( mtBank(jTube), gaussParams)
        end

        obj.featureBank = mtBank;

    end
    % }}}

    % EstimateMicrotubules {{{
    function obj = EstimateMicrotubules(obj, imXYZT, estParams)
        % EstimateMicrotubulePositions: estimates mt positions in a 2D/3D image.

        if obj.dim == 2               
            disp( 'Estimating microtubule locations...' )
        elseif obj.dim == 3
            disp( 'Estimating microtubule locations using a 2D MIP in the 3rd-dimension...' )
        end
        
        % extract the full cell image, and the segmented cell image
        imRegion = max( obj.sourceImage, [], 3);
        imCell = max( obj.sourceImage .* obj.mask, [], 3);
        imMask = logical( max( obj.mask, [], 3) );

        % Image Editing
        [imBP, imEst, imDisp] = EditMicrotubuleImageForEstimation( obj, imXYZT, imRegion, imMask);
        
        % Search for high intensity lines
        [coords, orients] = searchForHighIntensityLines3( imBP, imEst, imMask, estParams.interphase.estimation);
        obj.numberOfMicrotubules = length( coords);
        
        % Initialize microtubules with coordinate information
        obj = initializeMicrotubules( obj, coords, orients, imCell, imBP.*imMask);        
        obj.display.image = imDisp;

        % Estimate polynomial curve coefficients
        for jTube = 1 : obj.numberOfMicrotubules
            obj.featureBank( jTube) = EstimateMicrotubuleCurve( obj.featureBank( jTube), estParams.interphase.polyOrder);
        end
       
        PlotMicrotubuleBank( obj, 'estimatedCoords')
%         PlotMicrotubuleBank( obj, 'estimatedCoefs')
        return

        % Estimate gaussian parameters
        for jTube = 1 : numTubes
            obj.featureBank( jTube) = obj.featureBank( jTube).EstimateGaussianParameters( obj.featureBank(jTube) )
        end

    end
    % }}}

    % Fitting
    % fitMicrotubules {{{
    function obj = fitMicrotubules( obj, fitParams)

        % The microtubules are prepared for fitting
        % params( init, ub, lb) are passed, the problem is created, plotting function is defined, output function is created to add new parameters to file
        for jTube = 1 : obj.numberOfMicrotubules
            
            obj.featureBank( jTube) = Microtubule.fitInitialization( obj, fitParams.initialization)
            obj.featureBank( jTube) = Microtubule.fitMicrotubule( obj, fitParams.fitting)

        end

        % 

    end
    % }}}

end


% Secondary Methods
methods

    % [SWL, imGaussianMasked] = EditMicrotubuleImageForEstimation(imRegion) {{{
    function [imFilt, imSteer, imDisp] = EditMicrotubuleImageForEstimation(obj, imXYZT, imXY, imMask)

        disp('Enhancing the image...')

        % We will enhance the image differently for use in:
        % 1) Estimation/Detection process
        % 2) Display process

        %  Start by setting any pixels that are equal to 0 to the median value. This takes care of cells that have 0-padding on their edges.
        imXY( imXY == 0) = median( imXY(imXY ~= 0) );
        imXYZT( imXYZT == 0) = median( imXYZT(imXYZT ~= 0) );
        
        % Enhancement
        % We use a bandpass filter to only allow features of standard deviation between 1 and 4 (it removes noise variations and large scale changes in the cytoplasm background)
        sig_highpass = 1;
        sig_lowpass = 4;
        imFilt = filter_gauss_bandpass( imXY, sig_lowpass, sig_highpass);

        % Estimation Process
        % we take a mean in the z-dimension, this evens out the noise.  we run a steerable detector on each T-frame of imXYZT
        imXYT = squeeze( mean( imXYZT, 3) );
        imSteerZT = 0*imXYZT;
        imSteerT = 0*imXYT;

        for jT = 1 : size(imXYT, 3)
           
           for jZ = 1 : size(imXYZT, 3)


                imfilt{jZ, jT} = filter_gauss_bandpass( imXYZT(:,:,jZ,jT), sig_lowpass, sig_highpass);
                [~, ~, imSteerZT(:,:,jZ,jT)] = steerableDetector( imfilt{jZ,jT}, 4, sqrt(2)*1.3 );
            end
            imSteerT = squeeze( sum( imSteerZT, 3) );

        end

        imSteer = sum(imSteerT, 3).*imerode(imMask, strel('disk', 1) );
        

        T_steer = multithresh(imSteer, 1);
        imSteer( imSteer < T_steer(1) ) = 0;

        imSteer = imSteer.*imFilt;

        % Display
        imDisp = imFilt; imDisp( imDisp == 0) = median( imDisp( imDisp ~= 0) );
        T = multithresh(imDisp( imMask), 3); 
        imDisp( imDisp < T(1)) = 0;
        imDisp = imadjust( mat2gray(imDisp), [T(1), 1 ] );
        imDisp = mat2gray( imgaussfilt( imDisp, 1) ).*imMask;
        
        
        dispImg( imFilt.*imMask, imDisp, imSteer, [1 3]);
        title( sprintf('Cell %d: bandpass, display, steerableFilter' , obj.sourceInfo.cellNumber)) 
        set(gcf, 'Name', ['cell_' num2str(obj.sourceInfo.cellNumber) '_ImgEnh'], 'NumberTitle', 'off')
        
    end
    % }}}

    % imMicrotubulesCell = ParseMicrotubulesFromNetwork( imNetwork, imageIntensity, plotflag) {{{
    function imMicrotubulesCell = ParseMicrotubulesFromNetwork( obj, imNetwork, imageIntensity, imMask, guessParams)

        % Image Understanding 
        disp( 'Understanding the image...' )

        imBombed = bwareafilt( imNetwork, [10, Inf]); 
        imBombed = destroyNetworkJunctions( imBombed);
        imBombed = bwareafilt( imBombed, [10, Inf]);
        
        cc = bwconncomp( imBombed);
        for jBone = 1 : cc.NumObjects
                
            imBone = 0*imBombed;
            imBone( cc.PixelIdxList{jBone} ) = 1;
                
            skelProps(jBone) = analyzeLineImage( imBone);
                
        end

        if cc.NumObjects == 0
            imMicrotubulesCell = {};
            return
        end
        % cause breaks in the skeleton whenever there is a sharp orientation change
        % in an object.
        xRm = []; yRm = []; 
        for jj = 1 : length(skelProps)
            signal = smooth( skelProps(jj).thetaRaw, 3 );
        
            s1 = signal; s2 = signal;
            s1(end-1:end) = []; s2( 1:2)=[];
            x = skelProps(jj).orderedPixelList(:, 2); x(1) = []; x(end) = [];
            y = skelProps(jj).orderedPixelList(:, 1); y(1) = []; y(end) = []; 
            deriv3 = s1-s2;
            idxKill = find( abs(deriv3) > 0.3);
            xRm = [xRm; x(idxKill) ];
            yRm = [yRm; y(idxKill) ];
        end 

        idxKill = sub2ind( size(imBombed), yRm, xRm);
        imTemp = imBombed; 
        imTemp( idxKill) = 0;
        imBroken = bwareafilt( imTemp, [5, Inf]);

        stD = regionprops( imdilate(imBroken, strel('disk', 1) ), imageIntensity, 'MeanIntensity', 'PixelIdxList');
        st = regionprops( imBroken, imageIntensity, 'MeanIntensity', 'PixelIdxList');
        idxAvoid = [ cat( 1, stD(:).PixelIdxList); find(imMask == 0) ];
        vals = imageIntensity(:); vals( idxAvoid) = NaN;
        mu = median( vals( ~isnan(vals) ) );
        sig = std( vals( ~isnan(vals) ) );
        objKeep = find([st.MeanIntensity] > mu+sig);
        imNew = 0*imageIntensity;
        for jj = 1 : length(objKeep)
            imNew( st( objKeep(jj) ).PixelIdxList ) = 1;
        end

        im2 = repmat( mat2gray( imageIntensity).* ~imNetwork, 1, 1, 3); im2(:,:,3) = 0*imageIntensity; im2( :,:,2) = imNetwork;
        im3 = repmat( mat2gray( imageIntensity).* ~imNew, 1, 1, 3); im3(:,:,3) = 0*imageIntensity; im3( :,:,2) = imNew;        
        
        
        % Now lets reconstruct the network by connecting any obviously connected
        % pieces. We put a strick constraint on orientation matching but allow
        % linking to objects as far as 15 pixels (~1.5 micron)

        % define parameters for network reconstruction
        costs.maxCost = 3;
        costs.maxDistLink = 20;
        costs.maxDistLinkPerp = 2;
        costs.maxPhiDiff_EE = deg2rad( 20);
        costs.maxPhiDiff_EL = deg2rad(10);
        costs.EEvsELScalingFactor = 1.5;
        costs.linkObjectSizeScalingFactor = 1.2;
        costs.linkObjectSizeForScaling = 10;
        paramsReconstruct.costs = costs;
        clear costs
        costs.maxPhiDiff_EE = deg2rad( 70);
        costs.maxDistLink = 10;
        paramsFindComp.costs = costs;
        
        % Reconstruct MT network
        imConn = reconstructNetwork(imNew, imageIntensity, paramsReconstruct);
    
        im4 = repmat( mat2gray( imageIntensity).* ~imConn, 1, 1, 3); im4(:,:,3) = 0*imageIntensity; im4( :,:,2) = imConn;
        if guessParams.interphase.plotflag.estimation_networkCleaning
            dispImg( imageIntensity, im2, im3, im4, [ 2 2]); set(gcf, 'Name','cleaned_up_network', 'NumberTitle', 'off')
        end

        imHelper = findMinimalComponents( logical(imConn), paramsFindComp, guessParams.interphase.plotflag.estimation_mtIsolation);

        rmSmallLines = 1;

        if rmSmallLines
            fprintf('Number of estimated lines is %d. Removing Small Lines\n', length( imHelper) );
            minLineSize = 10;
            rmLines = cellfun( @(v) sum( v(:) )<minLineSize, imHelper);
            imHelper( rmLines) = [];
        end
        imMicrotubulesCell = imHelper;

    end
% }}}

% mtoc = findMTOC( logicalImageCell, intensityImage, plotflag) {{{
    function mtoc = findMTOC( obj, logicalImageCell, imDisplay, estParams)
        % intensity image could be a single image or a cell containing multiple images
        
        % stop if no MTOCs found.
        if length( logicalImageCell ) == 0
            mtoc = [];
            disp('No Tracks and no MTOCs were found! No microtubules can exist for this cell.')
        return
        end
        
        % It is also possible that the widest point along the line is an MTOC,
        % and it has an overlap of tubes going different ways. We'll look for the
        % brightest pixel in a convolved image.
        if ~iscell( logicalImageCell)
            logicalImageCell = {logicalImageCell};
        end
        
        im2D = max( obj.sourceImage, [], 3);
        imMask = max( obj.mask, [], 3);

        imGauss = mat2gray( imgaussfilt( im2D .* imMask, 2) );
        for jLine = 1 : length( logicalImageCell)
            
            
            cLine = logicalImageCell{ jLine };

            [ ~, mtoc( jLine ).idx ] = max(imGauss(:) .* cLine(:) );
  
            [ mtoc( jLine ).y , mtoc( jLine ).x ] = ind2sub( size(im2D), mtoc(jLine).idx ); % note y and x are reversed because first index refers to row number, i.e the vertical axis
    
            % Lets also find the direction of nucleation
            boneProps = analyzeLineImage( cLine);
            mtoc( jLine ).phi = boneProps.thetaRaw( boneProps.orderedIdxList == mtoc( jLine).idx );
    
        end

        % plot MTOCs on top of masked image
        if estParams.interphase.plotflag.estimation_mtoc 
            imRGB = mat2gray(repmat( imDisplay, 1, 1, 3)); imRGB( :,:, 2:3) = 0;
            dispImg( imRGB); hold on
            plot( [mtoc.x], [mtoc.y], 'w*', 'MarkerSize', 12, 'LineWidth', 3); hold off;
            title('MTOC locations')
            set(gcf, 'Name', 'mtoc_locations', 'NumberTitle', 'off')
        end

    end

% }}}

% mtBank = initializeMicrotubules( imMicrotubulesCell, imCell, mtoc) {{{
    function obj = initializeMicrotubules( obj, coordsMT, orients, imOriginal, imDisplay)
        
        colors = distinguishable_colors( obj.numberOfMicrotubules, {'w', 'k'} );
        for jMT = 1 : obj.numberOfMicrotubules
            
            mtBank( jMT) = Microtubule( imOriginal, imDisplay, coordsMT{jMT}, orients( jMT) );
            mtBank( jMT).id = jMT;
            mtBank( jMT).display.color = colors( jMT, :);

        end
        if obj.numberOfMicrotubules == 0
                mtBank = [];
        end
        obj.featureBank = mtBank;

    end
% }}}

% idxDead = removeDeadMicrotubules( mtBank) {{{
    function idxDead  = flagDeadMicrotubules( mtBank)

        idxDead = find( [mtBank.dead]);
        
    end
% }}}

end

% Display/Plotting Methods
methods
        % PlotMicrotubuleBank {{{
        function PlotMicrotubuleBank( obj, type)
        % PlotMicrotubuleBank: used for all kinds of bank plots ( i.e combining all Microtubules into a single plot)
        %  type = { 'estimatedCoords', 'estimatedCoefs', 'fitCoefs'}
            
            if ~strcmp( type, 'estimatedCoords') && ~strcmp( type, 'estimatedCoefs') && ~strcmp( type, 'fitCoefs')
                error( 'InterphaseMicrotubuleBank.PlotMicrotubuleBank : 2nd argument doesnt match possible strings allowed.' )
            end
            
            preT = [ 'Cell ' num2str(obj.sourceInfo.cellNumber) ': '];
            postT = ['_' num2str(obj.sourceInfo.cellNumber)]; 
            if strcmp( type, 'estimatedCoords')
                figName = sprintf( ['estimate_mtbank_points' postT]);
                figTitle = sprintf( [preT, 'Microtubule Point Estimates']);
            elseif strcmp( type, 'estimatedCoefs')
                figName = sprintf( ['estimate_mtbank_coefs' postT]);
                figTitle = sprintf( [preT, 'Microtubule Curve Estimates']);
            elseif strcmp( type, 'fitCoefs')
                figName = sprintf( ['fit_mtbank_coefs' postT]);
                figTitle = sprintf( [preT, 'Microtubule Curve Fits']);
            end

            figure('Name', figName, 'NumberTitle', 'off')
            pos = get(gcf, 'position');
            set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
            subplot(121); imagesc( max(obj.sourceImage.*obj.mask,[],3)); axis equal; colormap gray; title([preT, 'Original Image'])
            set(gca, 'xlim', [1 size( obj.display.image ,1)], 'ylim', [1 size( obj.display.image ,2)]); set(gca, 'FontSize', 15);
            subplot(122); imagesc( obj.display.image); axis equal; colormap gray;
            set(gca, 'xlim', [1 size( obj.display.image,1)], 'ylim', [1 size( obj.display.image,2)]); hold on;
            LH = []; L = {};

            for jTube = 1 : obj.numberOfMicrotubules 
                
                if strcmp( type, 'estimatedCoef')

                t = linspace(0,1);
                px = obj.featureBank( jTube).estimatedCoef(1,:);
                py = obj.featureBank( jTube).estimatedCoef(2,:);
                x = polyval( px, t );
                y = polyval( py, t );
                
                elseif strcmp( type, 'fitCoef')

                t = linspace(0,1);
                px = obj.featureBank( jTube).fitCoef(1,:);
                py = obj.featureBank( jTube).fitCoef(2,:);
                x = polyval( px, t );
                y = polyval( py, t );

                else

                x = obj.featureBank( jTube).estimatedCoords.x;
                y = obj.featureBank( jTube).estimatedCoords.y;
                
                end
                
                c = obj.featureBank( jTube).display.color;
                lw = obj.featureBank( jTube).display.LineWidth+2;
                ms = obj.featureBank( jTube).display.MarkerSize;
                id = obj.featureBank( jTube).id;

                % plot all estimated points with links in between them
                plot( x, y, '-', 'Color', c, 'LineWidth', lw, 'Marker', '.', 'MarkerSize', ms)

                % plot mtoc location
                plot( x(1), y(1), 'r*', 'LineWidth', lw/2, 'MarkerSize', ms+4)

                % create legend entry
                LH(jTube) = plot(nan, nan, '-', 'Color', c, 'MarkerSize', ms, 'LineWidth', lw); L{jTube} = ['MT ', num2str(id)];    

            end
            hold off
            legend(LH, L);
            set(gca, 'FontSize', 15);
            title( figTitle);

        end
        % }}}

end

% estimate locations of microtubules. this is the guess function

% add microtubules to this organizing center

% assign unique identification to each mi9crotubule

% fitting for each microtubuel

% global optmization of mnicrotubules

% tracks microtubules by linking detections over time

% handles all display properties


























end
