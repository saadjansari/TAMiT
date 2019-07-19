classdef SPB < Spot

    properties
        numMicrotubules 
    end

    methods ( Access = public )
       
        % SPB {{{
        function obj = SPB( position, amplitude, sigma, dim, ref_image)

            obj = obj@Spot( position, amplitude, sigma, dim, ref_image, 'SPB');

        end
        % }}}

    end

end
