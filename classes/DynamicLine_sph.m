classdef DynamicLine_sph < DynamicFeature
   properties
        position % (numPts x numDim x numTimes)
        startPosition % (numPts x numDim x numTimes)
        len
        theta
   end
    
   methods
       % Construct the object
       function obj = DynamicLine_sph( t0, t1, pos0, leng, theta, amp, color, time_step, size_voxels)
           
           obj = obj@DynamicFeature('Line', t0, t1, amp, color, time_step, size_voxels);
           %obj.position = pos;
           obj.startPosition = pos0;
           obj.len = leng;
           obj.theta = theta;
           
       end
       
       % Display feature
       function displayFeature( obj, ax, t_current)
           
           if ~obj.existAtTime(t_current)
               error('obj does not exist at this time')
           end
           pp = obj.position(:,:, t_current-obj.time_start+1 );
           line( pp(:,1), pp(:,2), 'color', obj.color, 'LineWidth', 2 )

       end
       
       % get Length
       function [lens, lens_err, times] = getLength(obj)
           
            times = obj.time_step*([obj.time_start:obj.time_end]);
            lens = zeros( size(times) );
            lens_err = lens;
            % For each time, get matched_feature, then get its coordinates,
            % and then compute length
            for jf = 1 : length( obj.matched_feats)
                cf = obj.matched_feats{jf};
                if isempty(cf)
                    lens(jf) = nan;
                    continue
                end
                
                dXYZ = obj.size_voxels .* ( cf.endPosition - cf.startPosition );
                errXYZdiff = sqrt( cf.err_startPosition.^2 + cf.err_endPosition.^2);
                errXYZ = obj.size_voxels .* errXYZdiff;
               
                lens(jf) = sqrt( sum( dXYZ.^2 ) );
                lens_err(jf) = (1/lens(jf) )*sqrt( sum( (dXYZ.*errXYZ).^2) );
            end
       end
       
       % get Amplitude
       function [amp, amp_err, times] = getAmplitude(obj)
           amp = obj.amplitude(1,:,1);
           amp_err = obj.amplitude(1,:,2);
           times = obj.time_step*([obj.time_start:obj.time_end]);
       end
       
       function obj = matchFeatures( obj, feats_det)
            % Match feature to detection features
            % Fills out unknown parameters in tracked feature
           
            obj.matched_feats = cell(1, +obj.time_end - obj.time_start);
            % Loop over times of existence
            for jt = 1 : 1+obj.time_end - obj.time_start
                
                % Get the possible detected features at this time
                possible_lines = feats_det{jt}.featureList{1}.featureList(2:end);
                
                % Get this feature's data at this time
                c_len = obj.len(1,jt,1);
                c_theta = obj.theta(:,jt,1)';
                
                % Loop over possible detected curves and find a match
                match_found = 0;
                eps = 1e-3;
                for jdc = 1 : length( possible_lines)
                    if abs(c_len - possible_lines{jdc}.length) < eps && ...
                            all( abs(c_theta - possible_lines{jdc}.theta) < eps )
                        obj.matched_feats{jt} = possible_lines{jdc};
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