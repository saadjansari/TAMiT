function [imNoisy, mainObj, params] = simCell(noise, type)

    addpath('classes');

    % Try to load parameter file
%     try
%         param = load('initSim.m');
%     end

    params.type = type;
    params.dim = 3;
    params.size = [150 150 7];
    params.ampMax = 1;
    params.sigMax = [2 2 1.5];
    params.sigMin = [1.2 1.2 1.0];
    params.sigMinSPB = [1.5 1.5 1.0];
    params.lenMax = 3;
    params.lenMin = 0.5;
    params.lenMu = 1;
    params.voxelSize = [0.1 0.1 0.5];
    params.noise = noise;
%     params.thetaMin = 0.3*pi;
%     params.thetaMax = 0.8*pi;

    % Simulate first frame
    switch params.type
        case 'Monopolar'

            % Get number of features between 1 and 5
            nLines = 4;

            % Get coordinates, amplitude, sigma of SPB
            c_spb = [75 75 4];
            a_spb = params.ampMax.*rand(1);
            sig_spb = getRand(params, params.sigMax, params.sigMinSPB);
            % Create spb
            Pole = Spot( c_spb, a_spb, sig_spb, params.dim, {}, {'Color', 'red'}); 
            Pole.label = 'spb';
             
            % Get coordinates, amplitude, sigma of lines 
            lines = {};
            
            % amplitude, sigma of lines - fixed
            a_mt = 0.3*a_spb + 0.5*a_spb.*rand(1);
            sig_mt = getRand(params, sig_spb, params.sigMin);
            for j = 1 : nLines

                % Get end coordinates, 
%                 len_mt = params.lenMin + (params.lenMax - params.lenMin).*rand(1);
%                 len_mt = exprnd(params.lenMu);
                len_mt = 3;
                phi_mt = 2*pi.*rand(1);
                theta_mt = pi*rand(1);
                
                % get end position, and reduce length until end point is
                % within image size
                outside = 1;
                while outside
                    if params.dim == 3
                        c_mt = c_spb + [ len_mt* [sin(theta_mt)*cos(phi_mt), sin(theta_mt)*sin(phi_mt), cos(theta_mt)]]./params.voxelSize;
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
                newMT = Line( c_spb, c_mt, a_mt, sig_mt, params.dim, {}, {'Color', 'blue', 'Linewidth', 2});
                lines = {lines{:}, newMT};
                
            end
            AsterObj  = AsterMT( params.dim, Pole, lines{:} );
            mainObj = MonopolarAster( params.dim, zeros(params.size), {AsterObj}, {});
            
        case 'Bipolar'
            error('not set up')

        case 'Interphase'
            error('not set up')

    end
    imgSim = mainObj.simulateFeature( params.size);
    imNoisy = addPoissNoise(imgSim, a_mt, params.noise);
%     figure; subplot(122); imagesc( max( imNoisy,[],3)); colormap gray; axis equal;
%     subplot(121); imagesc( max(imgSim,[],3)); colormap gray; axis equal;
%     figure; imagesc( max( imNoisy,[],3)); colormap gray; axis equal; xlim([0 150]); ylim([0 150]); xticks([]); yticks([]);
%     title(['Noise = ', num2str(params.noise)]);

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
    function imNoisy = addPoissNoise(img, ampMT, mu)
        noisyy = mu*ampMT*poissrnd( 50, size(img))/50;
        imNoisy = img + noisyy;
    end
    % }}}
end
