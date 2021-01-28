classdef DynamicFeatureMgmt
   properties
       features
       data_xyt
   end
    
   methods
       % Construct the object
       function obj = DynamicFeatureMgmt( features)
       
           obj.features = features;
           
       end
       
       % get lifetimes
       function lifetimes = getLifetimes( obj)
           lifetimes = zeros( 1, length(obj.features) );
           for jf = 1 : length(lifetimes)
               lifetimes(jf) = obj.features{jf}.getLifetime();
           end
       end
       
       % get lengths
       function [lens,lens_err,times] = getLengths( obj)
           lens = cell( 1, length(obj.features) );
           lens_err = cell( 1, length(obj.features) );
           times = cell( 1, length(obj.features) );
           for jf = 1 : length( lens)
               [lens{jf}, lens_err{jf}, times{jf}] = obj.features{jf}.getLength();
           end
       end
       
       % get amplitudes
       function [amp, amp_err, times] = getAmplitudes( obj)
           amp = cell( 1, length(obj.features) );
           amp_err = cell( 1, length(obj.features) );
           times = cell( 1, length(obj.features) );
           for jf = 1 : length( amp)
               [amp{jf}, amp_err{jf}, times{jf}] = obj.features{jf}.getAmplitude();
           end
       end
       
       % save data to file
       function saveMat(obj, path)
           tau = obj.getLifetimes();
           [lens,lens_err, times] = obj.getLengths();
           [amp,amp_err, ~] = obj.getAmplitudes();
           save([path, filesep,'dydata.mat'],'tau','lens', 'lens_err', 'amp', 'amp_err', 'times');
       end
      
       
   end
    
    
end