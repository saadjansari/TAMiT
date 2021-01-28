classdef DynamicFeature
   properties
       time_start
       time_end
       position % (numPts x numDim x numTimes)
       amplitude
       color
       time_step
       size_voxels
   end
    
   methods
       % Construct the object
       function obj = DynamicFeature( t0, t1, pos, amp, color, time_step, size_voxels)
       
           obj.time_start = t0;
           obj.time_end = t1;
           obj.position = pos;
           obj.amplitude = amp;
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
       function [lens, lens_err, times] = getLength(obj)
           lens = zeros( 1, obj.time_end-obj.time_start+1 );
           lens_err = lens;
           times = obj.time_step*([obj.time_start:obj.time_end]);
           for idx = 1 : length(lens)
               dXYZ = obj.size_voxels .* ( obj.position(2,:,idx,1) - obj.position(1,:,idx,1) );
               errXYZ = obj.size_voxels .* ( obj.position(2,:,idx,2) );
               
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
       
   end
    
    
end