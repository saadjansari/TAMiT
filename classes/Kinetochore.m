classdef Kinetochore < Spot

    methods ( Access = public )
       
        % Kinetochore {{{
        function obj = Kinetochore( position, amplitude, sigma, dim, ref_image)

            obj = obj@Spot( position, amplitude, sigma, dim, ref_image, 'Kinetochore');

        end
        % }}}
        
    end

end
