classdef DynamicFeature
    properties
        type
        time_start
        time_end
        amplitude
        color
        time_step
        size_voxels
        matched_feats
    end

    methods
       % Construct the object
       function obj = DynamicFeature( type, t0, t1, amp, color, time_step, size_voxels)
           
           if ~strcmp(type, 'Line') && ~strcmp( type, 'Curve') && ~strcmp( type, 'Spot')
                error('DynamicFeature: unknown type, allowed values are Line, Curve, Spot')
           end

           obj.type = type;
           obj.time_start = t0;
           obj.time_end = t1;
           obj.amplitude = amp;
           obj.color = color;
           obj.time_step = time_step;
           obj.size_voxels = size_voxels;

       end

       % Check if feature exists at a given time
       function status = existAtTime( obj, t_probe)
           if (t_probe >= obj.time_start) && (t_probe <= obj.time_end)
               status = true;
           else
               status = false;
           end
       end

       % get Lifetime
       function lifetime = getLifetime(obj)
           lifetime = obj.time_step * (obj.time_end - obj.time_start + 1);
       end

       % get Amplitude
       function [amp, amp_err, times] = getAmplitude(obj)
           amp = obj.amplitude(1,:,1);
           amp_err = obj.amplitude(1,:,2);
           times = obj.time_step*([obj.time_start:obj.time_end]);
       end

    end
    
    
end