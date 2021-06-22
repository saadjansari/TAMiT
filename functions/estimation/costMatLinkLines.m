function costMat = costMatLinkLines( movieInfo1, movieInfo2)

    % Function Layout:
        % Analyze Lines to find properties used for cost calculation
        % Define cost parameters for linking
        % Calculate cost of linking


    % Analyze the lines and obtain useful feature information {{{
    num1 = movieInfo1.num;
    num2 = movieInfo2.num;

    for j1 = 1 : num1 
        px = movieInfo1.polyCoef( j1, :, 1);
        py = movieInfo1.polyCoef( j1, :, 2);
        t = linspace( 0, 1);
        mts1(j1).com = [ polyval(px, 0.5); polyval(py, 0.5) ];
        mts1(j1).orient = atan( diff(polyval(py, [0,1]))/diff(polyval(px, [0,1])) );
        mts1(j1).endPt = [polyval( px, [0,1]) ; polyval(py, [0,1]) ];
        mts1(j1).length = norm( [px(1), py(1) ]);
        mts1(j1).coefX = px;
        mts1(j1).coefY = py;
    end
    for j1 = 1 : num2 
        px = movieInfo1.polyCoef( j1, :, 1);
        py = movieInfo1.polyCoef( j1, :, 2);
        t = linspace( 0, 1);
        mts2(j1).com = [ polyval(px, 0.5); polyval(py, 0.5) ];
        mts2(j1).orient = atan( diff(polyval(py, [0,1]))/diff(polyval(px, [0,1])) );
        mts2(j1).endPt = [polyval( px, [0,1]) ; polyval(py, [0,1]) ];
        mts2(j1).length = norm( [px(1), py(1) ]);
        mts2(j1).coefX = px;
        mts2(j1).coefY = py;
    end
    % }}}

    % Link Parameters {{{
    cost_scale = 5;
    axis_dist_max = 2;
    axis_phi_max = deg2rad( 15);
    coi_dist_max =  30;
    displayLinkInfo = 0;
    % }}}

    % Cost Calculation {{{
    cost_scale = 2;
    for j1 = 1 : num1, for j2 = 1: num2
    
        % 1. Axis Dist
        % find distance of mt1 from axis of mt2 and vice versa. axis dist is the mean of the two.
        axis_dist = findMinimumDistanceBetweenLines( mts1(j1), mts2(j2), 2 );
        costAxisDist = cost_axis_dist( axis_dist, axis_dist_max, cost_scale, 1);

        % 2. Axis Phi 
        % find difference in orientations of axis of the 2 microtubules
        axis_phi = abs( mts1(j1).orient - mts2(j2).orient);
        costAxisPhi = cost_axis_phi( axis_phi, axis_phi_max, cost_scale, 1);
       
        % 3. COI Dist
        % find distance from COI_1 to COI_2
        coi_dist = norm( mts1(j1).com - mts2(j2).com );
        costCoiDist = cost_coi_dist( coi_dist, coi_dist_max, 0.5*cost_scale, 1);
       
        costMat( j1, j2) = costAxisDist;

        if displayLinkInfo
            disp( sprintf( 'link: %d to %d', j1, j2) )
            disp( sprintf( 'axis_dist = %.2f', axis_dist) )
%             disp( sprintf( 'axis_phi= %.2f', axis_phi) )
%             disp( sprintf( 'coi_dist = %.2f', coi_dist) )
%             disp( sprintf( 'cost final = %.2f', costMat( j1, j2) ) )
            disp( '----------------------------')
        end

    end; end
    % }}}

    % findMinimumDistanceBetweenLines {{{
    function axisDist = findMinimumDistanceBetweenLines( mt1, mt2, measure)
        % We'll use the frechet or hausdorff distance between the curves to determine their similarity.
        % Using the discrete frechet or hausdorff measure

        t = linspace(0,1);
        % Get discrete points on curves
        x1 = polyval( mt1.coefX, t);
        y1 = polyval( mt1.coefY, t);
        x2 = polyval( mt2.coefX, t);
        y2 = polyval( mt2.coefY, t);
   
        if measure == 1
            % Frechet
            axisDist = DiscreteFrechetDist( [x1 ; y1], [x2 ; y2] );
        elseif measure == 2
            % Hausdorff
            axisDist = HausdorffDist( [x1' , y1'], [x2' , y2'] );
        end
        
        visualize = 0;
        if visualize == 1
        % temporary plot
        emptyfig = randn(150, 150);
        figure; imagesc( emptyfig); colormap gray; axis equal; hold on
        plot( x1, y1, 'r-'); plot(x2, y2, 'b-'); 
        hold off; title( sprintf('Axis Dist = %f', axisDist) )
        end

    end
    % }}}

    % Cost Calculation Functions: {{{
    % Measurements for Cost Calculation
    % 1. Axis Dist
    % 2. Axis Phi
    % 3. COI Dist

    function cost = cost_axis_dist(axis_dist, axis_dist_max, cost_scale, order)
        cost = ( cost_scale/ axis_dist_max^order ) * axis_dist^order;
    end

    function cost = cost_axis_phi(axis_phi, axis_phi_max, cost_scale, order)
        cost = ( cost_scale/ axis_phi_max^order) * axis_phi^order;
    end

    function cost = cost_coi_dist(coi_dist, coi_dist_max, cost_scale, order)
        cost = ( cost_scale/ coi_dist_max^order) * coi_dist^order;
    end
    % }}}

end
