function [imgSim, mainObj, params] = simMonopolar()

    addpath('../../classes');

    params.dim = 3;
%     params.snr = snr;
    params.sigMax = [2 2 1.5];
    params.sigMin = [1.2 1.2 1.0];
    params.sigMinSPB = [1.5 1.5 1.0];

    % Simulate first frame
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
    a_spb = 2 + 2*rand(1);
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
       
    imgSim = mainObj.simulateFeature( params.size);
%    imNoisy = addPoissNoise(imgSim, a_mt, snr);
%     figure; subplot(122); imagesc( max( imNoisy,[],3)); colormap gray; axis equal;
%     subplot(121); imagesc( max(imgSim,[],3)); colormap gray; axis equal;
%     figure; imagesc( max( imNoisy,[],3)); colormap gray; axis equal; xlim([0 params.size(1)]); ylim([0 params.size(1)]); xticks([]); yticks([]);
%     title(['SNR = ', num2str(params.snr)]);

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
    
end
