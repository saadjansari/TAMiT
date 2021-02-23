classdef DynamicLine < DynamicFeature
   properties
       position % (numPts x numDim x numTimes)
   end
    
   methods
       % Construct the object
       function obj = DynamicLine( t0, t1, pos, amp, color, time_step, size_voxels)
           
           obj = obj@DynamicFeature('Line', t0, t1, amp, color, time_step, size_voxels);
           obj.position = pos;
           
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
           lens = zeros( 1, obj.time_end-obj.time_start+1 );
           lens_err = lens;
           times = obj.time_step*([obj.time_start:obj.time_end]);
           for idx = 1 : length(lens)
               dXYZ = obj.size_voxels .* ( obj.position(2,:,idx,1) - obj.position(1,:,idx,1) );
               errXYZdiff = sqrt( obj.position(1,:,idx,2).^2 + obj.position(2,:,idx,2).^2);
               errXYZ = obj.size_voxels .* errXYZdiff;
               
               lens(idx) = sqrt( sum( dXYZ.^2 ) );
               lens_err(idx) = (1/lens(idx) )*sqrt( sum( (dXYZ.*errXYZ).^2) );
               
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
                c_endpos = obj.position(2,:,jt,1);
                
                % Loop over possible detected curves and find a match
                match_found = 0;
                eps = 1e-3;
                for jdc = 1 : length( possible_lines)
                    if all( abs(c_endpos - possible_lines{jdc}.endPosition) < eps)
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