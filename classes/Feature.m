classdef Feature < handle & matlab.mixin.Copyable
% This is a Feature superclass, it will be utilized for the Microtubule and the Kinetochore subclasses, with room for additional feature specializations.
    properties
%         position
%         amplitude
%         sigma
        dim
        image
        imageSim
%         display % display properties
        type % feature specialization (if any)
%         props2Fit
        ID 
    end

    methods

        % Feature {{{
        function obj = Feature( dim, image, type)
        % Feature : this is the constructor function for a Feature

            % Ensure dim matches image dimensionality
            if dim ~= length( size( image) )
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end

            obj.dim = dim;
            obj.image = image;
            
            % generic feature if unspecified, and speicalized feature is specified.
            if nargin < 3, 
                obj.type = 'generic';
            else 
                obj.type = type; 
            end

            % Assign a unique ID at random
            obj.ID = java.rmi.server.UID();

        end
        % }}}
        
        function objCopy = copyShallow(obj)
            % creates a shallow copy
            objCopy = copy( obj);
        end
    end

end
