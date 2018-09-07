classdef InterphaseMicrotubuleBank
% defines an organizing center for interphase microtubules. The bank can be thought of as a set of rules that governs the interactions between individual features. The bank is sort of an administrator at a given time, that figures out how individual features act when viewed together.
properties
    sourceImage
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
    function obj = InterphaseMicrotubuleBank( image3D, mask, movieName, imageTime, cellLocation, cellLocationImage, cellNumber)
    
        % All that is needed to initialize a bank is a source image( could be 2D or 3D).         
        dim = length( size(image3D));
        giveError( {dim, dim}, {'~=', '~='}, {2,3}, 'or', 'InterphaseMicrotubuleBank: 1st argument must be a 2D or a 3D image.')

        disp( 'Interphase Microtubule Bank created. ' )
        
        obj.sourceImage = image3D;
        obj.mask = mask
        obj.dim = length( size(image3D) );
        obj.sourceInfo.movieName = movieName;
        obj.sourceInfo.cellNumber = cellNumber;
        obj.sourceInfo.time = imageTime;
        obj.sourceInfo.cellLocationInMovie = cellLocation;
        obj.sourceInfo.cellLocationImage = cellLocationImage;

    end
    % }}}

    % Estimation
    % EstimateMicrotubules {{{
    function obj = EstimateMicrotubules(obj, gaussParams)
        % EstimateMicrotubulePositions: estimates mt positions in a 2D/3D image.

        if obj.dim == 2               
            disp( 'Estimating microtubule locations...' )
        elseif obj.dim == 3
            disp( 'Estimating microtubule locations using a 2D MIP in the 3rd-dimension...' )
        end
        
        % extract the full cell image, and the segmented cell image
        imRegion = max( obj.sourceImage, [], 3);
        imCell = max( obj.sourceImage .* obj.mask, [], 3);
        imMask = max( obj.mask, [], 3);

        % Image Editing
        [imNetwork, imageGaussian] = EditMicrotubuleImageForEstimation( imRegion);
       
        % Parse the image to separate out microtubules
        imMicrotubulesCell = ParseMicrotubulesFromNetwork( imNetwork, imageGaussian, obj.display.plotflags.estimation.microtubuleIsolation);
         
        % find microtubule nucleation point by looking for the brightest intensity point along the microtubule.
        mtoc = findMTOC( imMicrotubulesCell, imCell, obj.display.plotflags.estimation.mtoc);
        
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
    function [SWL, imGaussianMasked] = EditMicrotubuleImageForEstimation(imRegion)

        disp('Editing the image...')                                                                                                                                                     
        % sets any 0-padding for edge cells to the median value.
        imEstimate = imRegion; imEstimate( imEstimate == 0) = median( imEstimate( imEstimate ~= 0) );
        T = multithresh(imEstimate( maskForCell), 3); 
        imEstimate( imEstimate < T(1)) = 0;
        imEstimate = imadjust( mat2gray(imEstimate), [T(1), 1 ] );
        imGaussian = mat2gray( imgaussfilt( imEstimate, 1.5) );
        imGaussianMasked = mat2gray( imGaussian .* imMask);

        [~, ~, nonmax] = steerableDetector( imGaussian, 4, 2); SWL = nonmax .* imMask;
        T3 = multithresh( SWL( SWL > 0), 4); SWL( SWL < T3(1) ) = 0; SWL = logical(SWL);     

    end
    % }}}

    % imMicrotubulesCell = ParseMicrotubulesFromNetwork( imNetwork, imageIntensity, plotflag) {{{
    function imMicrotubulesCell = ParseMicrotubulesFromNetwork( imNetwork, imageIntensity, plotflag)

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
        idxAvoid = [ cat( 1, stD(:).PixelIdxList); find(maskForCell == 0) ];
        vals = imGaussian(:); vals( idxAvoid) = NaN;
        mu = median( vals( ~isnan(vals) ) );
        sig = std( vals( ~isnan(vals) ) );
        objKeep = find([st.MeanIntensity] > mu+sig);
        imNew = 0*imCell;
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
        paramsFindComp.costs = costs;

        % Reconstruct MT network
        imConn = reconstructNetwork(imNew, imageIntensity, paramsReconstruct);
    
        im4 = repmat( mat2gray( imageIntensity).* ~imConn, 1, 1, 3); im4(:,:,3) = 0*imageIntensity; im4( :,:,2) = imConn;
        if plotflag
            dispImg( imCell, im2, im3, im4, [ 2 2]); set(gcf, 'Name','cleaned_up_network', 'NumberTitle', 'off')
        end

        imHelper = findMinimalComponents( logical(imConn), paramsFindComp, plotflag);

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
    function mtoc = findMTOC( logicalImageCell, intensityImage, plotflag)
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

        imGauss = mat2gray( imgaussfilt( intensityImage, 2) );
        for jLine = 1 : length( logicalImageCell)
    
            cLine = logicalImageCell{ jLine };
            [ ~, mtoc( jLine ).idx ] = max(imGauss(:) .* cLine(:) );
  
            [ mtoc( jLine ).y , mtoc( jLine ).x ] = ind2sub( size(intensityImage), mtoc(jLine).idx ); % note y and x are reversed because first index refers to row number, i.e the vertical axis
    
            % Lets also find the direction of nucleation
            boneProps = analyzeLineImage( cLine);
            mtoc( jLine ).phi = boneProps.thetaRaw( boneProps.orderedIdxList == mtoc( jLine).idx );
    
        end

        % plot MTOCs on top of masked image
        if plotflag
            imRGB = mat2gray(repmat( imGauss, 1, 1, 3)); imRGB( :,:, 2:3) = 0;
            dispImg( imRGB); hold on
            plot( [mtoc.x], [mtoc.y], 'w*', 'MarkerSize', 12, 'LineWidth', 3); hold off;
            title('MTOC locations')
            set(gcf, 'Name', 'mtoc_locations', 'NumberTitle', 'off')
        end

    end

% }}}

% mtBank = initializeMicrotubules( imMicrotubulesCell, imCell, mtoc) {{{
    function mtBank = initializeMicrotubules( imMicrotubulesCell, imCell, mtoc)

        % Ensure 1 MTOC, 2 microtubules per track. Assign helper images and correct
        % phi to each of the 2 microtubules per track.
        for jNuc = 1 : length( mtoc)
            imt = 1 + 2*(jNuc-1);
            fmt = imt+1;

            for jTube = imt:fmt
                if mod( jTube, 2)
                    mtBank( jTube) = Microtubule.Microtubule( imCell, mtoc( jNuc).x, mtoc( jNuc).y, mtoc( jNuc).phi );
                else
                    mtBank( jTube) = Microtubule.Microtubule( imCell, mtoc( jNuc).x, mtoc( jNuc).y, mtoc( jNuc).phi - pi );
                end
                % Assign the helper image
                mtBank( jTube).helperImage = mat2gray( imgaussfilt( 1.0*imdilate( imMicrotubulesCell{ jNuc}, strel('disk', 1)) .* imCell ) );
            end
        end

        % Assign identification and color to each microtubule. 
        numTubes = length( mtBank);
        colors = distinguishable_colors( numTubes, {'w', 'k'} );
        for jTube = 1 : numTubes
            mtBank( jTube).id = jTube;
            mtBank( jTube).display.color = colors( jTube, :);
        end

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

    



end

% estimate locations of microtubules. this is the guess function

% add microtubules to this organizing center

% assign unique identification to each mi9crotubule

% fitting for each microtubuel

% global optmization of mnicrotubules

% tracks microtubules by linking detections over time

% handles all display properties


























end
