function spindles = find_spindle( frames)
    
    % Params
    p.spindleDeterminationSensitivity = 0.4;
    p.spindleMinIntensity = 0.65;
    p.linewidth = 2;
    p.brightestPixelAsSPB = 0;
    
    %[f, ha] = disp_frames( 2, length(frames)/2, frames);
    
    masks = {};
    spindles = {};
    for jframe = 1 : length( frames)
        
        % Get frame
        frame = frames{jframe};
        mask = find_spindle_mask( frame, p);        
        masks = {masks{:}, mask};
        
        % find spindle endpoints
        spindle = find_spindle_line( frame, mask, p );
        
        % construct spindle object
        lineAmpList = Cell.findAmplitudeAlongLine( frame, spindle.MT.startPosition, spindle.MT.endPosition );
        spindleAmp = median(lineAmpList) - median( frame( frame(:) ~= 0) );
        spindleMT = Line( spindle.MT.startPosition, spindle.MT.endPosition, ...
            spindleAmp, [2,2,1.2], 3, {}, {'Color', [1,0,1], 'LineWidth',2});
        spindleMT.label = 'spindle'; spindleMT.SetBounds();
                
        spindles = {spindles{:}, spindleMT};
        
    end

    function mask = find_spindle_mask( frame, p)
        
        % Convolve the image with a 3D gaussian to bring out the signal from the SPB
        image3DConv = imgaussfilt3( frame, 1, 'FilterDomain', 'spatial');
        imPlane = mat2gray( max( image3DConv, [], 3) );

        % Keep the strongest signal pixels
        threshOtsu = max( [thresholdOtsu( imPlane( imPlane > 0) ) 0.4]);
        imPlaneStrong = imPlane;
        imPlaneStrong( imPlaneStrong < threshOtsu) = 0;

        %imMask = imextendedmax( imPlaneStrong, p.spindleDeterminationSensitivity);
        imPlaneStrongBool = imPlaneStrong;
        imPlaneStrongBool( imPlaneStrongBool > 0 ) = 1;
        mask = imPlaneStrongBool;
        
    end

    function spindle = find_spindle_line( frame, mask, p, ax)
        
        imPlane = max(frame,[],3);
        numVoxelsY = size( imPlane, 1); numVoxelsX = size( imPlane, 2);
        
        % we'll pick the maxima which has the biggest mean intensity
        SS = regionprops( mask, imPlane, 'Centroid', 'MajorAxisLength', 'Orientation', 'MeanIntensity');
        SScell = struct2cell( SS);
        meanInts = cell2mat( SScell(4,:) );
        [ spindleInt, idx] = max( meanInts);

        if spindleInt < 0.5
            error( 'The mean spindle intensity is less than 50% of the maximum intensity in the image');
        else
            fprintf('\t Success: Spindle Mean Intensity = %.2f\n', spindleInt)
        end

        % Now we can obtain the shape parameters of the spindle
        centroid = cell2mat( SScell( 1, idx) );
        majAxisLen = cell2mat(  SScell( 2, idx) );
        angleOrient = deg2rad( cell2mat(  SScell( 3, idx) ) );

        if angleOrient < 0
            angle = 3*pi/2 - abs(angleOrient);
        elseif angleOrient >=0
            angle = pi/2 + angleOrient;
        end
        
        % Find coordinates of a long line passing through the region
        if majAxisLen < 3
            majAxisLen = 3;
        end

        ptMin = round( flip(centroid) + (5*majAxisLen) * [ cos(angle), sin(angle)  ] );
        ptMax = round( flip(centroid) - (5*majAxisLen) * [ cos(angle), sin(angle)  ] );

        coords1 = round( linspace( ptMin(1), ptMax(1), ceil(12*majAxisLen) ) );
        coords2 = round( linspace( ptMin(2), ptMax(2), ceil(12*majAxisLen) ) );
        coords = [coords1 ; coords2];

        % remove any column in coords whose either entry is less than or equal
        % to 0. or greater than or equal to numVoxelsX or numVoxelsY
        rmCol = [];
        for jCol = 1 : length(coords)

            if any( find( coords( :, jCol) < 1) )
                rmCol = [ rmCol, jCol];
            elseif coords(1, jCol) > numVoxelsY
                rmCol = [ rmCol, jCol];
            elseif coords(2, jCol) > numVoxelsX
                rmCol = [ rmCol, jCol];
            end

        end
        % remove bad columns
        coords( :, rmCol) = [];
        
        % Find Intensity on line with a non-zero width
        perpMatrix = [-p.linewidth : p.linewidth] .* [cos( angle+pi/2 ) ; sin( angle+pi/2)];

        % measure intensity
        IntSpindle = zeros(1, length(coords) );
        idxMaxPix = zeros(1, length(coords) );
        imPlaneG = imgaussfilt( imPlane, 1);
        for jPix = 1 : length(IntSpindle)
            coordsPerp = round( coords(:, jPix) +  perpMatrix);
            coordsPerp( 1, find( coordsPerp(1,:) > numVoxelsY) ) = numVoxelsY; coordsPerp( 1, find( coordsPerp(1,:) < 1) ) = 1;
            coordsPerp( 2, find( coordsPerp(2,:) > numVoxelsX) ) = numVoxelsX; coordsPerp( 2, find( coordsPerp(2,:) < 1) ) = 1;
            idxPerp = sub2ind( [numVoxelsY, numVoxelsX], coordsPerp(1, :), coordsPerp(2, :) );
            [IntSpindle( jPix), idxMaxInt] = max( imPlaneG(idxPerp) );
            idxMaxPix( jPix) = idxPerp( idxMaxInt);
        end
        IntSpindle = mat2gray( smooth( IntSpindle, 4) );
        
        % Now lets obtain the two SPB points. 
        indSpindle = find( IntSpindle > p.spindleMinIntensity);
        if length( indSpindle) == 1
            error('spindle found is not long enough. spindle min intensity is possibly too high. ');
        end
        ind1 = indSpindle(1); ind2 = indSpindle(end);

        if length(indSpindle) > 7
            ind1 = ind1+1;
            ind2 = ind2-1;
        end

        % Here X is the horizontal coordinate.
        [AstersY, AstersX] = ind2sub( size(imPlane), idxMaxPix( [ind1, ind2] ) );

        if AstersX(1) == AstersX(2) && AstersY(1) == AstersY(2)
            error('only one point found for spindle')
        end
        
        if p.brightestPixelAsSPB

            % Now find the maximum intensity location in this image plane and set
            % this as one of the two SPBs
            [ ~, indMax] = max( imPlane(:) );
            [SPB1y, SPB1x] = ind2sub( size(imPlane), indMax);

            % We'll keep this max intensity SPB, and for the second we'll choose
            % one from the two we found by picking the one furthest from this max
            % position.
            dist1 = norm( [ SPB1y, SPB1x] - [ AstersY(1), AstersX(1)] );
            dist2 = norm( [ SPB1y, SPB1x] - [ AstersY(2), AstersX(2)] );

            if dist1 > dist2
                SPB2y = AstersY(1);
                SPB2x = AstersX(1);
            else
                SPB2y = AstersY(2);
                SPB2x = AstersX(2);
            end

            AstersX = [SPB1x, SPB2x];
            AstersY = [SPB1y, SPB2y];

        end
        
         % find Z coordinates of the spindle
        [~, zIdx] = max( imgaussfilt3( frame,1), [], 3);
        AstersZ = [ zIdx( round(AstersY(1)), round(AstersX(1)) ), zIdx( round(AstersY(2)), round(AstersX(2)) ) ];
        spindlePosition = [AstersX' , AstersY', AstersZ'];

        spindle.MT.startPosition = spindlePosition(1, :);
        spindle.MT.endPosition = spindlePosition(2, :);
        spindle.MT.amplitude = mean( Cell.findAmplitudeAlongLine( frame, spindlePosition(1, :), spindlePosition(2, :) ) );
        
        if nargin == 4
            axes( ax)
            line( AstersX, AstersY, 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 8, 'Color', 'r' )
        end
        
    end
    
end