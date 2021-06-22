function SpindleObjs = find_astrals( frames, spindles)
    
    [f, ha] = disp_frames( 2, length(frames)/2, frames);
    
    % Params
    p.spindleExclusionRange = deg2rad(45);
        
    SpindleObjs = {};
    for jframe = 1 : length( frames)
        
        % Get frame, spindle
        frame = frames{jframe};
        spindleMT = spindles{jframe};
        
        % Find the angle of this spindle w.r.t to each pole
        spindleAngle(1) = mod( atan2( spindleMT.endPosition(2)-spindleMT.startPosition(2) , spindleMT.endPosition(1)-spindleMT.startPosition(1) ) , 2*pi );
        spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

        % Get SPBs
        AsterObjects = create_spb_objects( frame, spindleMT);
        
        % Get start position of astrals (3-5 pixels away from the
        % spb opposite to the direction of spindle.
        rr = 2;
        sp{1} = spindleMT.startPosition + rr*[cos(spindleAngle(1)+pi), sin(spindleAngle(1)+pi), 0];
        sp{2} = spindleMT.endPosition + rr*[cos(spindleAngle(2)+pi), sin(spindleAngle(2)+pi), 0];
        
        % find phi Values
        phiVals1 = find_phi_peaks( frame, sp{1}, spindleAngle(1), p);
        phiVals2 = find_phi_peaks( frame, sp{2}, spindleAngle(2), p);

        cm{1} = find_lines_given_phi( frame, sp{1}, phiVals1, spindleAngle(1), p);
        cm{2} = find_lines_given_phi( frame, sp{2}, phiVals2, spindleAngle(2), p);
        
        for jAster = 1:2
            if ~isempty( cm{jAster} )
                AsterObjects{jAster} = { AsterObjects{jAster}{:}, cm{jAster}{:} };
            end
        end
        Asters = cell(1,2);
        for jAster = 1 : 2
            Asters{jAster} = Aster( 3, AsterObjects{jAster}{:} );
        end
        SpindleObj = SpindleNew( 3, frame, {spindleMT, Asters{:} }, {});
        SpindleObjs = { SpindleObjs{:}, SpindleObj };
        
        % Display
        axes( ha(jframe) ); hold on;
        SpindleObj.displayFeature( ha(jframe) ); 
        
    end

    function peakPhi = find_phi_peaks( imageIn, startPoint, spindleAngle, p, ax)
        
        % Do a sweep radial sweep and look for peaks corresponding to
        % rods.
        im2 = max( mat2gray(imageIn),[],3); 
        [~,~,nms,~] = steerableDetector( im2, 4, 2);
        im2g = imgaussfilt(im2, 1);
        [phiIntensity, phiValues] = Cell.radIntegrate2D( im2g, startPoint, 5, 15);

        % find peaks in Phi Intensity
        imVals = im2g( im2g ~= 0);
        minPkHeight = 1.0*median( imVals );

        warning('off', 'signal:findpeaks:largeMinPeakHeight' )
        [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, ...
            'MinPeakHeight', minPkHeight, 'MinPeakProminence', 0.25*median( imVals ));
        peakPhi = mod( peakPhi, 2*pi);
        warning('on', 'signal:findpeaks:largeMinPeakHeight' )
            
        % remove angle belonging to spindle.
        spindleAngle = mod( spindleAngle, 2*pi);
        idxRm = find( abs(peakPhi-spindleAngle) < p.spindleExclusionRange | ...
            abs(peakPhi-spindleAngle+2*pi) < p.spindleExclusionRange | ...
            abs(peakPhi-spindleAngle-2*pi) < p.spindleExclusionRange);
        peakPhi( idxRm)=[];

        % plot lines in phi directions
        if nargin == 5
            axes( ax)
            for jPhi = peakPhi
                plot( startPoint(1), startPoint(2), 'b*', 'MarkerSize', 6)
                plot( startPoint(1) + [0,10]*cos(jPhi), startPoint(2) + [0,10]*sin(jPhi), 'g-', 'LineWidth',2);
            end
        end
        
    end

    function mask = create_astral_mask( frame)
       
        imageIn = mat2gray( frame);
        
        % create a mask to apply to steerable image
        mask = imgaussfilt3(imageIn, 2);
        st = Methods.GetImageStats(mask,0);
        % mask( mask < st.Median+2*st.Sigma) = 0; mask(mask ~=0) = 1;
        mask( mask < 1.0*st.Median) = 0; mask(mask ~=0) = 1;
        
    end

    function curvedMTs = find_lines_given_phi( frame, startPoint, peakPhi, spindleAngle, p, ax)
        
        % Default Values
        defaultVisibility = 10;
        defaultFieldOfView = 70;
        defaultStepSize = 4;
        minLength = 10;
        
        imageIn = mat2gray( frame);
        nms3 = 0*imageIn;
        for jZ = 1 : size(imageIn,3)
            [~, ~, nms3(:,:,jZ), ~] = steerableDetector( imageIn(:,:,jZ), 4, 3);
        end
        mask = create_astral_mask( frame);
        AstralMicrotubules = {};
            
        for pp = 1 : length(peakPhi)

            imSearch = sum(nms3,3).*max(mask,[],3);
            % Find the max intensity start point in close neighborhood
            sx = startPoint(1); sy = startPoint(2); spp = [sx,sy];
            maxInt = imSearch( round(sy), round(sx));
            for jx = sx+[-1,0,1]
                for jy = sy+[-1,0,1]
                    int = imSearch( round(jy), round(jx));
                    if int > maxInt
                        maxInt = int;
                        spp = [jx,jy];
                    end
                end
            end

            % Set the two opposite angles for search
            cc = Methods.estimateCurveCoords( spp', peakPhi(pp), imSearch, defaultStepSize, defaultVisibility, defaultFieldOfView, 1);
            cc(1:2,1) = startPoint(1:2)';
            
            % figure out z-coordinates
            [a2D, i2d] = max( imageIn, [], 3);
            c3 = zeros( 1, size(cc,2));
            for jj = 1 : length(c3)
                rng = 1; alist = []; idlist = [];
                for j1 = -rng : rng
                    for j2 = -rng : rng

                        if round(cc(2,jj))+j1 > 1 && round(cc(2,jj))+j1 < size(mask,1) && ...
                                round(cc(1,jj))+j2 > 1 && round(cc(1,jj))+j2 < size(mask,2)
                            alist = [ alist, a2D( round(cc(2,jj))+j1, round(cc(1,jj))+j2)];
                            idlist = [ idlist, i2d( round(cc(2,jj))+j1, round(cc(1,jj))+j2)];
                        end

                    end
                end
                [~, idd] = max( alist); c3( jj) = idlist(idd);
            end
                
            c3(1: ceil(length(c3)/2) ) = c3(1);

            try; c3(end) = mean( c3(end-2:end-1));end
                
            c3(ceil(length(c3)/2) : end) = mean( c3(ceil(2*length(c3)/3):end) );
            pz = polyfit( linspace(0,1, length(c3)), c3, 1);
            zz = polyval( pz, linspace(0,1, length(c3)));
            zz(zz < 1) = 1.2; zz(zz > size(imageIn,3) ) = size(imageIn,3)-0.2;
            cc = [cc; zz];

            % Length of end-end curve
            if norm( cc(:,end)-cc(:,1)) > minLength
                AstralMicrotubules = { AstralMicrotubules{:}, cc};
            end

        end
            
        % Remove any duplicates (also remove if angle matches the
        % spindle angle)
        phiz = zeros(1, length(AstralMicrotubules));
        for j1 = 1 : length( AstralMicrotubules)
            phiz(j1) = atan2( AstralMicrotubules{j1}(2,2)-AstralMicrotubules{j1}(2,1) , AstralMicrotubules{j1}(1,2)-AstralMicrotubules{j1}(1,1) );
        end

        % spindle angle first
        idxRm = find( abs(phiz-spindleAngle) < p.spindleExclusionRange | ...
            abs(phiz-spindleAngle+2*pi) < p.spindleExclusionRange | ...
            abs(phiz-spindleAngle-2*pi) < p.spindleExclusionRange);
        AstralMicrotubules(idxRm) = []; phiz(idxRm) = [];

        destroy = [];
        for p1 = 1:length(phiz)
            for p2 = 1:length(phiz)
                if p1 ~=p2
                    if abs( phiz(p1) - phiz(p2)) < 0.2 || abs( phiz(p1) - phiz(p2) -2*pi) < 0.2 || abs( phiz(p1) - phiz(p2) + 2*pi) < 0.2
                        % figure out which one to destroy
                        if norm( AstralMicrotubules{p1}(:,end)-AstralMicrotubules{p1}(:,1)) > norm( AstralMicrotubules{p2}(:,end)-AstralMicrotubules{p2}(:,1))
                            destroy = [destroy, p2];
                        else
                            destroy = [destroy, p1];
                        end
                    end
                end
            end
        end
        destroy = unique(destroy);
        AstralMicrotubules( destroy) = [];
        
        % create curved mts
        curvedMTs = create_curved_objects( frame, AstralMicrotubules);
        
        if ~isempty(curvedMTs)
            % ensure mts are within the image region
            curvedMTs = restrict_curvedMT_inside_image( frame, curvedMTs);

            % reduce mt length based on image intensity
            curvedMTs = threhold_mt_length_using_image( frame, curvedMTs);
        end
    end

    function curvedMTs = restrict_curvedMT_inside_image( frame, curvedMTs)
        
        nvx = size(frame,2);
        nvy = size(frame,1);
        nvz = size(frame,3);
        % restrict mts within the image region
        for jf = 1 : length(curvedMTs)
            
            outside = 0;
            % GetCoords
            cc = curvedMTs{jf}.GetCoords();
            
            % Identify coords that exceed the size of the image, and remove them
            rmIdx_x = find( cc(1,:) < 1 | cc(1,:) > nvx);
            rmIdx_y = find( cc(2,:) < 1 | cc(2,:) > nvy);
            rmIdx_z = find( cc(3,:) < 1 | cc(3,:) > nvz);
            rmIdx = union( union( rmIdx_x , rmIdx_y ), rmIdx_z);
            if ~isempty(rmIdx)
                outside = 1;
                disp('This line estimate exceeds the image region. Restricting it inside')
            end
            cnt = 0;
            while outside && cnt < 20
                cnt = cnt+1;
                curvedMTs{jf}.L = 0.95*curvedMTs{jf}.L;
                cc = curvedMTs{jf}.GetCoords();
                rmIdx_x = find( cc(1,:) < 1 | cc(1,:) > nvx);
                rmIdx_y = find( cc(2,:) < 1 | cc(2,:) > nvy);
                rmIdx_z = find( cc(3,:) < 1 | cc(3,:) > nvz);
                outside = ~isempty( union( union( rmIdx_x , rmIdx_y ), rmIdx_z) );
            end
            
        end
        
    end

    function curvedMTs = threhold_mt_length_using_image( frame, curvedMTs)
        % restrict mts within the image region
        
        frame = mat2gray(frame);
        %mask = create_astral_mask( frame);
        img = imgaussfilt(frame,1);
        %med = median( img( :) );
        
        threshold_factor = 1.1;
        
        %f = figure;
        %ha = tight_subplot( 1, length(curvedMTs), 0.01, 0.01,0.01);
        
        % For each curved mt, gets its coords
        for jf = 1 : length(curvedMTs)
            
            % GetCoords
            cc = curvedMTs{jf}.GetCoords();
            
            % Get intensity at these coordiantes
            int = Cell.findAmplitudeAlongCurveCoords( img, round( cc) );
            tt = 1 : length(int);
            
            % Find a mask for this specific feature by usign its simulated
            % image, then binarizing it, dilating it.
            mask = curvedMTs{jf}.simulateFeature( size(frame) );
            mask = imbinarize( mask, 0.005);
            mask = imdilate( mask, strel('disk',8,6) );
            %mask3 = repmat( mask, [1,1, size(frame,3)]);
            med = median( img( find( mask(:) ) ) );
            threshold = threshold_factor*med;
            
            idxEnd = length(int);
            if any( int < threshold)
                idxEnd = 1+length(int) - find( flip(int)> threshold, 1, 'first');
            end
            curvedMTs{jf}.L = (idxEnd/length(int) ) *curvedMTs{jf}.L;
            % plot things
%             axes( ha(jf) );
%             plot( tt, int, 'b-', 'LineWidth', 3); hold on;
%             plot( [tt(1) tt(end)], [med med], 'r:', 'LineWidth', 1.5)
%             plot( [tt(1) tt(end)], 1.5*[med med], 'm:', 'LineWidth', 1.5)
%             plot( [tt(1) tt(end)], 2*[med med], 'g:', 'LineWidth', 1.5)
        end
        
    end
    
    function curvedMTs = create_curved_objects( imageIn, AMT )
        
        imageIn = mat2gray(imageIn);
        
        bkg = median(imageIn( imageIn(:)~=0)); 
        
        idxRm = [];
        curvedMTs = cell(1, length( AMT));
        for jb = 1 : length( curvedMTs )

            coords = AMT{jb};
            % Get length
            L = sum( sqrt( diff( coords(1,:)).^2 + diff( coords(2,:)).^2 + diff( coords(3,:)).^2 ) );
            nInt1 = round( L/ length( coords(1,:) ) );
            [cX1,cY1,cZ1] = Methods.InterpolateCoords3( coords(1,:), coords(2,:), coords(3,:), nInt1 );
            % Get Coeff
            cf1 = Bundle.estimatePolyCoefficients( [cX1;cY1;cZ1], [3 3 1], linspace(0,L,length(cX1 )));
            % Get coordinates from coeffs
            t1 = linspace(0,L,length(cX1 ));
            x1 = polyval( cf1{1}, t1); y1 = polyval( cf1{2}, t1);

            % Get origin
            origin = coords(:,1);
            if origin(3) >= size( imageIn,3)
                origin(3) = size(imageIn,3)-0.2;
            elseif origin(3) <= 1
                origin(3) = 1.2;
            end
            % Get initial tangent vector and theta vector
            tanInit{1} = [cf1{1}(3), cf1{2}(3), cf1{3}(1)];
            thetaInit = [atan2( tanInit{1}(2), tanInit{1}(1) ), pi/2];
            % Normal Magnitude Coefficients
            nV = [ 2*(cf1{1}(3)*cf1{2}(2) - cf1{1}(2)*cf1{2}(3)), ...
                6*(cf1{1}(3)*cf1{2}(1) - cf1{1}(1)*cf1{2}(3))];

            % Get amplitude along each bundle
            A1 = smooth( Cell.findAmplitudeAlongCurveCoords( max(imageIn,[],3), round([cX1;cY1]) ) - bkg);
            A1( A1 < bkg) = bkg; amp = median(A1);

            % Ensure coefficients are within the image region
            thetaInit(1,2) = pi/2 - 0.03;
            if sign(nV(1)) == sign(nV(2)) && (abs(nV(1)) > 0.01 && abs(nV(2)) > 0.005)
                idxRm = [idxRm; jb];
            end

            % Create
            curvedMTs{jb} = CurvedMT( origin, thetaInit, nV, L, amp, [2,2,1.2], 3, {}, {'Color', 'red', 'LineWidth',2});
        end
        curvedMTs( idxRm) = [];        
        
    end

    function AsterObjects = create_spb_objects( frame, spindleMT)
        
        bkg = median( frame(frame(:)~=0) );
        
        % Spindle Poles
        AsterObjects{1} = {};
        AsterObjects{2} = {};
        
        % Create the SpindlePoleBody objects
        spbAmp(1) = frame( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) ) - bkg - spindleMT.amplitude;
        spbAmp(2) = frame( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) )-bkg-spindleMT.amplitude;
        if any(spbAmp < 0)
            spbAmp( spbAmp < 0) = bkg;
            warning( 'findFeaturesMT_deNovo : forcing SPBAmp to be > 0')
        end
        SPB{1} = Spot( spindleMT.startPosition, spbAmp(1), [3 3 1.2], 3, {}, {'Color', 'green', 'marker', 'o'}); 
        SPB{1}.label = 'spb'; SPB{1}.SetBounds();
        SPB{2} = Spot( spindleMT.endPosition, spbAmp(2), [3 3 1.2], 3, {}, {'Color', 'green', 'marker', 'o'}); 
        SPB{2}.label = 'spb'; SPB{2}.SetBounds();
        AsterObjects{1} = {SPB{1}};
        AsterObjects{2} = {SPB{2}};
        
    end

end