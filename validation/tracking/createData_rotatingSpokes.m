function createData_rotatingSpokes()

    % Initialize variables {{{
    dim = 2;
    origin = [0,0];
    numLines= 5;
    lineLength = 30;
    polyOrder = 1;
    omega = pi/500;
    totalT = 100;
    % }}}

    % Initialize Spokes {{{
    thetas0 = linspace( 0, 2*pi, numLines+1); thetas0(end) = [];
    for jt = 1 : numLines 
        % find final point by trig evaluations
        final = origin + [ cos( thetas0(jt) ) , sin( thetas0(jt) ) ];
        % fit polynomials to x and y
        px = polyfit( [0, 1], [final(1), origin(1)], polyOrder); 
        py = polyfit( [0, 1], [final(2), origin(2)], polyOrder); 
        poly(jt, 1).x = px; poly(jt, 1).y = py;
    end
    % }}}

    % Time Evolution and updating poly-coef {{{
    timeList = 1: totalT;
    delTheta = omega .* timeList;
    for jtime = timeList
        thetas = thetas0 + delTheta( jtime);
        for jt = 1 : numLines
            final = origin + [ cos( thetas(jt) ) , sin( thetas(jt) ) ];
            % fit polynomials to x and y
            px = polyfit( [0, 1], [final(1), origin(1)], polyOrder); 
            py = polyfit( [0, 1], [final(2), origin(2)], polyOrder); 
            poly( jt, jtime+1).x = px; poly( jt, jtime+1).y = py;
        end
    end
    % }}}

    % Plot spinning wheel
    %  initial frame
    figure; hold on
    for jt = 1 : numLines
        xx = polyval( poly( jt, 1).x, [0, 1] );
        yy = polyval( poly( jt, 1).y, [0, 1] );
        plot( xx, yy)
    end
    hold off; title( sprintf('Time = %f', 1) )
    drawnow
    pause(0.1)

    % future frames
    for jtime = timeList
        for jt = 1 : numLines
            xx = polyval( poly( jt, jtime+1).x, [0, 1] );
            yy = polyval( poly( jt, jtime+1).y, [0, 1] );
            plot( xx, yy); hold on;
        end
        hold off; title( sprintf('Time = %f', jtime) )
        drawnow
        pause(0.1)
    end

end
