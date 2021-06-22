
function newCurveModel()


    % initial values
    x0 = [0 0 0]';
    v0 = [-1 -1 0]';
    dt = 0.1;
    tmax = 2;
    t = 0:dt:tmax;

    % Acceleration : a(t) = a0 + a1*t + a2*t^2 
    % Case 1: Varying a0;
    a0 = linspace(-0.2,0.2,5); a1 = 0;
    
    figure; hold on;
    legs = {};
    for val = a0
        xx = GetCoords( x0, v0, t, dt, val, a1);
        plot( xx(1,:), xx(2,:), 'linewidth',3 );
        legs = {legs{:}, [ 'a(t) = ', num2str(val), '+0t']};
    end
    title('Varying a0')
    legend( legs{:})
    
    % Case 2: Varying a1;
    a0 = 0; a1 = linspace(-0.2,0.2,5);

    figure; hold on;
    legs = {};
    for val = a1
        xx = GetCoords( x0, v0, t, dt, a0, val);
        plot( xx(1,:), xx(2,:), 'linewidth',3 );
        legs = {legs{:}, [ 'a(t) = ', num2str(a0), '+('  num2str(val) , ')t']};
    end
    title('Varying a1')
    legend(legs{:})
    
    % Case 3: Varying t and a0;
    a0 = linspace(-0.2,0.2,5); a1 = 0;
    tmax = linspace(2,4,5);
    
    figure; hold on;
    legs = {};
    for jj = 1:length(a0)
        xx = GetCoords( x0, v0, 0:0.1:tmax(jj), dt, a0(jj), a1);
        plot( xx(1,:), xx(2,:) , 'linewidth',3);
        legs = {legs{:}, [ 'a0 = ', num2str(a0(jj)), ', t_{max} = ' num2str(tmax(jj))]};
    end
    title('Varying a0 and t')
    legend(legs{:})

    function xx = GetCoords( x0, v0, t, dt, a0, a1)
        
        % Acceleration function discretized in time
        acc = a0 + a1*t;

        % initialization
        xx = zeros( length( x0), length(t) ); xx(:,1) = x0;
        vv = zeros( 2, length(t) ); vv(:,1) = v0(1:2);
        Rot = [ 0 -1; 1 0]; % Rotation matrix

        % iterate
        for jt = 2 : length(t)
            vv(:,jt) = vv(:,jt-1) + acc( jt)*dt*Rot* vv(:,jt-1) / norm( vv(:,jt-1));
            v_mean = 1/2 * ( vv(:,jt) + vv(:,jt-1) );
            xx(1:2,jt) = xx(1:2,jt-1) + dt*v_mean/norm(v_mean);
            if length(x0)==3
                xx(3,jt) = xx(3,jt-1) + v0(3)*dt;
            end
        end
        
    end

end
    