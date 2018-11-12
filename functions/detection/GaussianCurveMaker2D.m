function imFeatures = GaussianCurveMaker2D( p, ImageFit, FixedParams)

%             save('tempCurveData.mat', 'p', 'ImageFit', 'FixedParams')
            speedUp = 2.0;
            nummt = FixedParams.numberOfMicrotubules; 
            fitdim = FixedParams.dim;
            order = FixedParams.PolyOrder;
            if fitdim == 3, orderZ = FixedParams.PolyOrderZ; end
            structFit = convertVec2Struct_MT( p, FixedParams);
            
            Bkg = structFit.background;
           
            % Initialize the image volume
            imPlane = 0*ImageFit;
            numPixX = size(imPlane, 1);
            numPixY = size(imPlane, 2);
            if fitdim == 3, numPixZ = size(imPlane, 3); end

            % load the meshgrids
            if fitdim == 2, [ yGrid, xGrid ] = meshgrid( 1 : numPixY , 1 : numPixX); xV = xGrid(:); yV = yGrid(:); gridCell={xV, yV};
            elseif fitdim == 3, [ yGrid, xGrid, zGrid ] = meshgrid( 1 : numPixY , 1 : numPixX, 1 : numPixZ ); xV = xGrid(:); yV = yGrid(:); zV = zGrid(:); gridCell = {xV, yV, zV}; end

            idx = 1 : numel( imPlane);

            % create parametric coord
            Tvec = linspace(0,1, round(100/speedUp) );
            Tvec( end) = [];

            % Now for each value of the parameter, we'll do a gauss
            % quadrature numerical integration 
            DeltaT = median( diff(Tvec) );

            % get the offsets for gauss quadrature
            poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
            poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;
            Tsort = sort([  Tvec + poff1,  Tvec + poff2] );
            
            imBkg = (imPlane + Bkg) .* FixedParams.Mask;
            imFeatures = imBkg;
            for jmt = 1 : nummt
                
                Amp = structFit.amplitude(jmt);
                Sig = structFit.std{jmt};

                XCoef = structFit.coefX{jmt};
                YCoef = structFit.coefY{jmt};
                % for speed.
                xLoc = polyval( XCoef, Tsort )';
                yLoc = polyval( YCoef, Tsort )';
                if fitdim==3, ZCoef = structFit.coefZ{jmt}; zLoc = polyval( ZCoef, Tsort)'; end
                
                dParam = 1;
                for jk = 1 : length( xLoc) -1
                    if fitdim == 2, dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2)];
                    elseif fitdim == 3, dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2 + diff(zLoc(jk:jk+1)).^2)]; end
                end
                dParamNew = linspace( dParam(1), dParam(end), length(dParam) );

                fitX = fit( dParam', xLoc, 'linearinterp');
                fitY = fit( dParam', yLoc, 'linearinterp');
                xDat = fitX( dParamNew);
                yDat = fitY( dParamNew);
                if fitdim==3, fitZ = fit( dParam', zLoc, 'linearinterp'); zDat = fitZ( dParamNew); end

                % if microtubule pokes out in z dimension out of real volume, then we set values to infinity.
                if fitdim==3 && ( any(zDat < 1) || any( zDat > numPixZ) )
                    imFeatures = inf( size(imPlane) );
                    warning('Z poked out of the real volume. Forced Infinities to guide lsqnonlin.')
                    return
                end

                if fitdim == 2, PtGaussLoc = [ xDat , yDat ];
                elseif fitdim == 3, PtGaussLoc = [ xDat , yDat , zDat]; end
                
    %             figure; imagesc(ImageFit); colormap gray, axis equal; hold on
    %             plot( xDat, yDat, 'r-', 'linewidth', 3); hold off

                % Compute the approximate integral over the spline points by
                % multiplying gaussians in x, y and z and adding for each
                % different value of the parameter.
                imConv = NumericalConv( Sig, PtGaussLoc, gridCell, imPlane, idx );
                imFeature = ( Amp - Bkg) * mat2gray( imConv);
                % set very small values to 0
                imFeature( imFeature < 0.1*max( imFeatures(:) ) ) = 0;
                imFeatures = imFeatures + imFeature;

            end

            % Penalize lsqnonlin if it draws a curve outside the mask region
            imCheck = imFeatures .* imcomplement( FixedParams.Mask);
            if any( imCheck(:) > 0.9*max(imFeatures(:) ) )
                imFeatures = imFeatures + Inf.*imcomplement( FixedParams.Mask);
            end
%             imPlane = imPlane .* FixedParams.Mask;           

            % NumericalConv {{{
            function imConv = NumericalConv( Sig, PtGaussLoc, gridIndexCell, imPlane, idx )

                % Create copies of image volume
                imConv = imPlane .* 0;

                fitdim = length(Sig); if length(gridIndexCell) ~= fitdim, error('Issue with passing grid indexes to NumericalConv2.m.'), end

                NumGauss = size( PtGaussLoc, 1);
                
                CoordsX = gridIndexCell{2} .* ones( 1, NumGauss);
                CoordsY = gridIndexCell{1} .* ones( 1, NumGauss);
                xx = PtGaussLoc( :, 1);
                yy = PtGaussLoc( :, 2);
                s1 = Sig(1);
                s2 = Sig(2);
                if fitdim==3,  
                    CoordsZ = gridIndexCell{3} .* ones( 1, NumGauss);
                    zz = PtGaussLoc( :, 1);
                    s3 = Sig(1);
                end

                Xvals = exp( -( xx' - CoordsX ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
                Yvals = exp( -( yy' - CoordsY ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);
                if fitdim == 2,
                    imConv( idx') = sum( Xvals .* Yvals , 2);
                elseif fitdim == 3
                    Zvals = exp( -( zz' - CoordsZ ).^2 / (2 * s3^2) ) / ( 2 * s3^2 * pi)^(1/2);
                    imConv( idx') = sum( Xvals .* Yvals .* Zvals , 2);
                end

            end
            % }}}
                    
end
