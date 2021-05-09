classdef DynamicSpot < DynamicFeature
    properties
        position
    end
    
    methods
        % Construct the object
        function obj = DynamicSpot( t0, t1, pos0, amp, color, time_step, size_voxels)

            obj = obj@DynamicFeature('Spot', t0, t1, amp, color, time_step, size_voxels);
            obj.position = pos0;

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

       
        function obj = matchFeatures( obj, feats_det)
            % Match feature to detection features
            % Fills out unknown parameters in tracked feature
           
            obj.matched_feats = cell(1, +obj.time_end - obj.time_start);
            % Loop over times of existence
            for jt = 1 : 1+obj.time_end - obj.time_start
                
                % Get the possible detected features at this time
                possible_spots = feats_det{jt}.featureList(1:end);
                
                % Get this feature's data at this time
                c_pos = obj.position(:,jt,1)';
                
                % Loop over possible detected curves and find a match
                match_found = 0;
                eps = 1e-3;
                for jdc = 1 : length( possible_spots)
                    if all( abs(c_pos - possible_spots{jdc}.position) < eps) 
                        obj.matched_feats{jt} = possible_spots{jdc};
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