function imPlane = GaussianCurveMaker2D( p, ImageFit, FixedParams)

%             save('tempCurveData.mat', 'p', 'ImageFit', 'FixedParams')
            nummt = FixedParams.numberOfMicrotubules; 
            fitdim = FixedParams.dim;
            order = FixedParams.PolyOrder;
            Bkg = p(1);

            structFit = convertVec2Struct( p, FixedParams);

            if FixedParams.StartingPointFixed == 1
                XCoef = [ p(5: 5+order-1), FixedParams.XCoefEnd ];
                YCoef = [ p(5+order: end), FixedParams.YCoefEnd ];
            elseif FixedParams.StartingPointFixed == 0
                XCoef = p(5: 5+order-1);
                YCoef = p(5+order: end);
            end
           
            % Initialize the image volume
            imPlane = 0*ImageFit;
            numPixX = size(imPlane, 1);
            numPixY = size(imPlane, 2);

            % load the meshgrids
            [ yGrid, xGrid ] = meshgrid( 1 : numPixY , 1 : numPixX);
            xV = xGrid(:); yV = yGrid(:);
            idx = 1 : numel( imPlane);

            % create parametric coord
            Tvec = linspace(0,1);
            Tvec( end) = [];

            % Now for each value of the parameter, we'll do a gauss
            % quadrature numerical integration 
            DeltaT = median( diff(Tvec) );

            % get the offsets for gauss quadrature
            poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
            poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;
            Tsort = sort([  Tvec + poff1,  Tvec + poff2] );
            
            for jmt = 1 : nummt
                
                XCoef = structFit.coefX{jmt};
                YCoef = structFit.coefY{jmt};
                if 
                % for speed.
                xLoc = polyval( XCoef, Tsort )';
                yLoc = polyval( YCoef, Tsort )';
                
                dParam = 1;
                for jk = 1 : length( xLoc) -1
                    dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2)];
                end

                fitX = fit( dParam', xLoc, 'linearinterp');
                fitY = fit( dParam', yLoc, 'linearinterp');
                dParamNew = linspace( dParam(1), dParam(end), length(dParam) );
                xDat = fitX( dParamNew);
                yDat = fitY( dParamNew);

                PtGaussLoc = [ xDat , yDat ];
                
    %             figure; imagesc(ImageFit); colormap gray, axis equal; hold on
    %             plot( xDat, yDat, 'r-', 'linewidth', 3); hold off

                % Compute the approximate integral over the spline points by
                % multiplying gaussians in x, y and z and adding for each
                % different value of the parameter.
                Conv2 = NumericalConv2( Sig, PtGaussLoc, xV, yV, imPlane, idx );
                
                imFeatures = ( Amp - Bkg) * mat2gray(Conv2);
                % set very small values to 0
                imFeatures ( imFeatures < 0.1*max( imFeatures(:) ) ) = 0;

            imBkg = (imPlane + Bkg) .* FixedParams.Mask;
            imPlane = imBkg + ( Amp - Bkg) * mat2gray( imFeatures);

            % Penalize lsqnonlin if it draws a curve outside the mask region
            imCheck = imPlane .* imcomplement( FixedParams.Mask);
            if any( imCheck(:) > 0.9*max(imFeatures(:) ) )
                imPlane = imPlane + Inf.*imcomplement( FixedParams.Mask);
            end
            imPlane = imPlane .* FixedParams.Mask;           
            
%             dispImg(ImageFit, imPlane, [1 2])
            % NumericalConv2 {{{
            function imConv2 = NumericalConv2( Sig, PtGaussLoc, xV, yV, imPlane, idx )

                % Create copies of image volume
                imConv2 = imPlane .* 0;

                NumGauss = size( PtGaussLoc, 1);

                CoordsX = yV .* ones( 1, NumGauss);
                CoordsY = xV .* ones( 1, NumGauss);

                xx = PtGaussLoc( :, 1);
                yy = PtGaussLoc( :, 2);

                s1 = Sig;
                s2 = Sig;

                Xvals = exp( -( xx' - CoordsX ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
                Yvals = exp( -( yy' - CoordsY ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);

                imConv2( idx') = sum( Xvals .* Yvals , 2);

            end
            % }}}
                    
end
