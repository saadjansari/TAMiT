classdef IMT_Manager < Feature

    properties
        nArrays
    end

    methods
        
        % IMT_Manager {{{
        function obj = IMT_Manager( position, amplitude, sigma, dim, ref_image)

            obj.position = position;
            obj.amplitude = amplitude;
            obj.sigma = sigma;
            obj.dim = dim;
            obj.ref_image = ref_image;
            
            obj = obj@Feature( position, amplitude, sigma, dim, ref_image, 'IMT_Manager');

        end
        % }}}

    end
end
