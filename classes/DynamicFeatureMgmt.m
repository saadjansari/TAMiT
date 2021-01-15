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
       function [lens,times] = getLengths( obj)
           lens = cell( 1, length(obj.features) );
           times = cell( 1, length(obj.features) );
           for jf = 1 : length( lens)
               [lens{jf}, times{jf}] = obj.features{jf}.getLength();
           end
       end
       
       % save data to file
       function saveMat(obj, path)
           tau = obj.getLifetimes();
           [lens,times] = obj.getLengths();
           save([path, filesep,'dydata.mat'],'tau','lens', 'times');
       end
      
       
   end
    
    
end