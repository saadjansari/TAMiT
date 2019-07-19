classdef MT_curved < Feature

    properties
        
    end

    methods
       
        % MT_curved {{{
        function obj = MT_Array( position, amplitude, sigma, dim, ref_image)

            obj.position = position;
            obj.amplitude = amplitude;
            obj.sigma = sigma;
            obj.dim = dim;
            obj.ref_image = ref_image;
            
            obj = obj@Feature( position, amplitude, sigma, dim, ref_image, 'MT_curved');

        end
        % }}}
        
    end
end
