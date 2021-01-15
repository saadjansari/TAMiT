classdef DynamicFeature
   properties
       time_start
       time_end
       position % (numPts x numDim x numTimes)
       color
       time_step
       size_voxels
   end
    
   methods
       % Construct the object
       function obj = DynamicFeature( t0, t1, pos, color, time_step, size_voxels)
       
           obj.time_start = t0;
           obj.time_end = t1;
           obj.position = pos;
           obj.color = color;
           obj.time_step = time_step;
           obj.size_voxels = size_voxels;
           
       end
       
       % Check if feature exists at a given time
       function status = existAtTime( obj, t_probe)
           if (t_probe >= time_start) && (t_probe <= time_end)
               status = true;
           else
               status = false;
           end
       end
       
       % Display feature
       function displayFeature( obj, ax, t_current)
           
           if ~obj.existAtTime(t_current)
               error('obj does not exist at this time')
           end
           pp = obj.position(:,:, t_current-obj.time_start+1 );
           line( pp(:,1), pp(:,2), 'color', obj.color, 'LineWidth', 2 )

       end
       
       % get Lifetime
       function lifetime = getLifetime(obj)
           lifetime = obj.time_step * (obj.time_end - obj.time_start + 1);
       end
       
       % get Length
       function [lens,times] = getLength(obj)
           lens = zeros( 1, obj.time_end-obj.time_start+1 );
           times = obj.time_step*([obj.time_start:obj.time_end]);
           for idx = 1 : length(lens)
               realCC = obj.size_voxels .* ( obj.position(2,:,idx) - obj.position(1,:,idx) );
               lens(idx) = sqrt( sum( realCC.^2 ) );
           end
       end
       
   end
    
    
end