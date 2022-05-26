classdef Feature < handle & matlab.mixin.Copyable
% This is a Feature superclass, it will be utilized for the Microtubule and the Kinetochore subclasses, with room for additional feature specializations.
    properties
        dim
        imageSim
        type % feature specialization (if any)
        ID
        label = ''
        fit = 'all' % 'all' or 'zonly'
    end

    methods

        % Feature {{{
        function obj = Feature( dim, type)
        % Feature : this is the constructor function for a Feature

            obj.dim = dim;
            
            % generic feature if unspecified, and speicalized feature is specified.
            if nargin < 2
                obj.type = 'generic';
            else 
                obj.type = type; 
            end

            % Assign a unique ID at random (Can be set by an organizer)
%             obj.ID = java.rmi.server.UID();

        end
        % }}}
        
        function objCopy = copyShallow(obj)
            % creates a shallow copy
            objCopy = copy( obj);
        end
        
        
    end

end
