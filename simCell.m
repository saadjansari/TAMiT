function [imNoisy, mainObj, params] = simCell(snr, type)

    addpath('classes');

    % Try to load parameter file
%     try
%         param = load('initSim.m');
%     end

    params.type = type;
    params.dim = 3;
    params.snr = snr;
    params.sigMax = [2 2 1.5];
    params.sigMin = [1.2 1.2 1.0];
    params.sigMinSPB = [1.5 1.5 1.0];
%     params.size = [120 120 15];



    % Simulate first frame
    switch params.type
        case 'Monopolar'
            %params
            params.voxelSize = [0.1 0.1 0.5];
            params.size = [80 80 7];
            zScale = params.voxelSize(3)/params.voxelSize(1);
            params.lenMax = 30;
            params.lenMin = 7;
            params.lenMu = 1;
            params.thetaMin = 0.3*pi;
            params.thetaMax = 0.8*pi;
    
            % Simulate pole/spb
            % coordinates, amplitude, sigma of SPB
            c_spb = [40 40 4];
            a_spb = 1 + 2*rand(1);
            sig_spb = getRand(params, params.sigMax, params.sigMinSPB);
            % create spb
            Pole = Spot( c_spb, a_spb, sig_spb, params.dim, {}, {'Color', 'red'}); 
            Pole.label = 'spb';
            
            % Simulate lines
            nLines = 4;
            lines = {};
            phiList = [];
            % amplitude, sigma of lines - fixed
            a_mt = 1;
            
            for j = 1 : nLines

                len_mt = (params.lenMax-params.lenMin)*rand(1) + params.lenMin;
                
                again=1;
                while again
                    phi_mt = 2*pi.*rand(1); 
                    if any( abs(phi_mt-phiList) < 0.5)
                        again=1;
                    else
                        phiList = [phiList, phi_mt];
                        again=0;
                    end
                end
                theta_mt = mod( acos( (-2*rand(1)+1)/zScale), pi);
%                 theta_mt = pi*rand(1);
                sig_mt = getRand(params, sig_spb, params.sigMin);
                % get end position, and reduce length until end point is
                % within image size
                outside = 1;
                while outside
                    if params.dim == 3
                        c_mt = c_spb + [ len_mt* [sin(theta_mt)*cos(phi_mt), sin(theta_mt)*sin(phi_mt), cos(theta_mt)]];
                    elseif params.dim == 2
                        c_mt = c_spb + len_mt* [cos(phi_mt), sin(phi_mt)];
                    end
                    if any(c_mt < 1) || any( c_mt > params.size)
                        outside = 1;
                        if 0.95*len_mt > params.lenMin
                            len_mt = 0.95*len_mt;
                        else
                            theta_mt = pi*rand(1);
                        end
                    else
                        outside = 0;
                    end
                end
                
                % Create line 
                disp(len_mt)
                newMT = Line( c_spb, c_mt, a_mt, sig_mt, params.dim, {}, {'Color', 'blue', 'Linewidth', 2});
                lines = {lines{:}, newMT};
                
            end
            AsterObj  = AsterMT( params.dim, Pole, lines{:} );
            mainObj = MonopolarAster( params.dim, zeros(params.size), {AsterObj}, {});
            
        case 'Bipolar'
            error('not set up')

        case 'Interphase'
            error('not set up')
            
        case 'MitosisBud'
            %params
            params.voxelSize = [0.05 0.05 0.5];
            params.size = [120 120 15];
            zScale = params.voxelSize(3)/params.voxelSize(1);
            params.lenMax = 60;
            params.lenMin = 7;
            params.lenMu = 1;
            params.thetaMin = 0.3*pi;
            params.thetaMax = 0.8*pi;
            
            % Simulate spindle
            % 2 poles + connecting line
            % Get length of spindle between 5 and 10 pixels
            lenSpindle = (10-5)*rand(1) + 5;
            thetaSpindle = pi/2;
            phiSpindle = 2*pi*rand(1);
            c1 = ceil(params.size/2) + (lenSpindle/2)*[sin(thetaSpindle)*cos(phiSpindle), ...
                sin(thetaSpindle)*sin(phiSpindle), cos(thetaSpindle)/zScale];
            c2 = ceil(params.size/2) - (lenSpindle/2)*[sin(thetaSpindle)*cos(phiSpindle), ...
                sin(thetaSpindle)*sin(phiSpindle), cos(thetaSpindle)/zScale];
            % coordinates, amplitude, sigma of SPB
            a_spb = 2*rand(1);
            sig_spb = getRand(params, params.sigMax, params.sigMinSPB);
            % create spb
            Pole{1} = Spot( c1, a_spb, sig_spb, params.dim, {}, {'Color', 'red'}); 
            Pole{1}.label = 'spb';
            Pole{2} = Spot( c2, a_spb, sig_spb, params.dim, {}, {'Color', 'red'}); 
            Pole{2}.label = 'spb';
            % create spindle line
            a_s = 3 + 2*rand(1);
            sig_s = getRand(params, params.sigMax, params.sigMinSPB);
            spindleMT = Line( c1, c2, a_s, sig_s, params.dim, {}, {'Color', 'red', 'LineWidth',4});
            spindleMT.label = 'spindle';
            
            spindleAngle(1) = mod( atan2( spindleMT.endPosition(2)-spindleMT.startPosition(2) , spindleMT.endPosition(1)-spindleMT.startPosition(1) ) , 2*pi );
            spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

            % Get start position of astrals (3-5 pixels away from the
            % spb opposite to the direction of spindle.
            rr = 2.5;
            sp{1} = spindleMT.startPosition + rr*[cos(spindleAngle(1)+pi), sin(spindleAngle(1)+pi), 0];
            sp{2} = spindleMT.endPosition + rr*[cos(spindleAngle(2)+pi), sin(spindleAngle(2)+pi), 0];
            
            Asters = cell(1,2);
            for jA = 1:2
                
                % Simulate lines
                nLines = 2;
                lines = {};
                
                % amplitude, sigma of lines - fixed
                a_mt = 1;
                for j = 1 : nLines

                    len_mt = (params.lenMax-params.lenMin)*rand(1) + params.lenMin;
                    thetaInit = zeros(1,2);
                    thetaInit(1) = mod( rand(1)+pi+spindleAngle(jA)-0.5, 2*pi);
                    thetaInit(2) = mod( acosd( ((cos(params.thetaMax)-cos(params.thetaMax))*rand(1) + cos(params.thetaMax))/zScale), 2*pi);
%                     thetaInit = mod( [ rand(1)+pi+spindleAngle(jA)-0.5 , (params.thetaMax - params.thetaMin)*rand(1) + params.thetaMin],2*pi);
                    if randi([0 1]) == 0
                        nV = [0.01*rand(1), -0.005*rand(1)];
                    else
                        nV = [-0.01*rand(1), 0.005*rand(1)];
                    end
                    sig_mt = getRand(params, sig_spb, params.sigMin);
                    % get params, and reduce length until end point is
                    % within image size
                    newMT = CurvedMT( sp{jA}, thetaInit, nV, len_mt, a_mt, sig_mt, params.dim, {}, {'Color', 'blue', 'Linewidth', 2});
                    cc = newMT.GetCoords();
                    while any( cc(:,end) < [1,1,1]') || any( cc(:,end) > params.size')
                        %reduce length until hits lower bound
                        if 0.95*len_mt > params.lenMin
                            len_mt = 0.95*len_mt;
                        else
                            thetaInit = [ thetaInit(1) , (params.thetaMax - params.thetaMin)*rand(1) + params.thetaMin];
                        end
                        newMT = CurvedMT( sp{jA}, thetaInit, nV, len_mt, a_mt, sig_mt, params.dim, {}, {'Color', 'blue', 'Linewidth', 2});
                        cc = newMT.GetCoords();
                    end
                    lines = {lines{:}, newMT};
                end
                Asters{jA}  = AsterMT( params.dim, Pole{jA}, lines{:} );
            end
            mainObj = SpindleNew( params.dim, zeros(params.size), {spindleMT, Asters{:}}, {});

    end
    imgSim = mainObj.simulateFeature( params.size);
    imNoisy = mat2gray(addPoissNoise(imgSim, a_mt, snr));
%     figure; subplot(122); imagesc( max( imNoisy,[],3)); colormap gray; axis equal;
%     subplot(121); imagesc( max(imgSim,[],3)); colormap gray; axis equal;
    figure; imagesc( max( imNoisy,[],3)); colormap gray; axis equal; xlim([0 params.size(1)]); ylim([0 params.size(1)]); xticks([]); yticks([]);
    title(['SNR = ', num2str(params.snr)]);

    % getCoordsRand {{{
    function cc = getRand(params, ub, lb)

        if nargin < 3
            lb = zeros( size(ub));
        end

        switch params.dim
            case 2
                % Get x, y, z
                x = ( ub(1) - lb(1) ).*rand(1) + lb(1);
                y = ( ub(2) - lb(2) ).*rand(1) + lb(2);
                cc = [x, y];
            case 3
                % Get x, y, z
                x = ( ub(1) - lb(1) ).*rand(1) + lb(1);
                y = ( ub(2) - lb(2) ).*rand(1) + lb(2);
                z = ( ub(3) - lb(3) ).*rand(1) + lb(3);
                cc = [x, y, z];
            otherwise 
                error('dim should be 2 or 3')
        end
    end
    % }}}

    % getCoords {{{
    function cc = getCoords(params)
        switch params.dim
            case 2
                % Get x, y, z
                x = params.size(1).*rand(1);
                y = params.size(2).*rand(1);
                cc = [x, y];
            case 3
                % Get x, y, z
                x = params.size(1).*rand(1);
                y = params.size(2).*rand(1);
                z = params.size(3).*rand(1);
                cc = [x, y, z];
            otherwise 
                error('dim should be 2 or 3')
        end
    end
    % }}}
    
    % addPoissNoise {{{
    function imNoisy = addPoissNoise(img, ampMT, snr)
        noisyy = poissrnd( 5, size(img))/5;
        imNoisy = img + (ampMT/snr)*noisyy;
    end
    % }}}
end
