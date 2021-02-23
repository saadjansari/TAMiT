classdef DynamicCurve < DynamicFeature
    properties
        startPosition % (numPts x numDim x numTimes)
        len
        theta
    end
    
    methods
        % Construct the object
        function obj = DynamicCurve( t0, t1, pos0, leng, theta, amp, color, time_step, size_voxels)

            obj = obj@DynamicFeature('Curve', t0, t1, amp, color, time_step, size_voxels);
            obj.startPosition = pos0;
            obj.len = leng;
            obj.theta = theta;

        end
       
        % Display feature
        function displayFeature( obj, ax, t_current)

            if ~obj.existAtTime(t_current)
               error('obj does not exist at this time')
            end
            error('not set up')
            %            pp = obj.position(:,:, t_current-obj.time_start+1 );
            %            line( pp(:,1), pp(:,2), 'color', obj.color, 'LineWidth', 2 )

        end
       
        % get Length
        function [lens, lens_err, times] = getLength(obj)
            
            times = obj.time_step*([obj.time_start:obj.time_end]);
            lens = zeros( size(times) );            
            % For each time, get matched_feature, then get its coordinates,
            % and then compute length
            for jf = 1 : length( obj.matched_feats)
                cf = obj.matched_feats{jf};
                if isempty(cf)
                    lens(jf) = nan;
                    continue
                end
                cc = cf.GetCoords();
                lens(jf) = sum( sqrt( sum( diff(cc.*obj.size_voxels',1, 2).^2,1) ) );
            end
           
            lens_err = zeros( size(obj.len) );
        end
       
        function obj = matchFeatures( obj, feats_det)
            % Match feature to detection features
            % Fills out unknown parameters in tracked feature
           
            obj.matched_feats = cell(1, +obj.time_end - obj.time_start);
            % Loop over times of existence
            for jt = 1 : 1+obj.time_end - obj.time_start
                
                % Get the possible detected features at this time
                possible_curves = feats_det{jt}.featureList(2:end);
                
                % Get this feature's data at this time
                c_len = obj.len(1,jt,1);
                c_theta = obj.theta(:,jt,1)';
                
                % Loop over possible detected curves and find a match
                match_found = 0;
                eps = 1e-3;
                for jdc = 1 : length( possible_curves)
                    if abs(c_len - possible_curves{jdc}.GetLength()) < eps && ...
                            all( abs(c_theta - possible_curves{jdc}.thetaInit) < eps )
                        obj.matched_feats{jt} = possible_curves{jdc};
                        obj.matched_feats{jt}.display{2} = obj.color;
                        match_found = 1;
                        break
                    end
                end
                if ~match_found
                    disp(['MATCH NOT FOUND - ',num2str(jt)])
                end
                    
            end
       end
       
   end
    
    
end