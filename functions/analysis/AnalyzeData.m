classdef AnalyzeData
    % This is the data analysis class. It is basically a container for all activities analysis related
    properties
        routine
        channel = []
        time = []
        path
    end

    methods (Access = public )

        % AnalyzeData {{{
        function obj = AnalyzeData( path, routine, varargin)
           
            if strcmp( routine, 'frame')
                obj.routine = routine;
            elseif strcmp( routine, 'channel')
                obj.routine = routine;
                if length( varargin) ~= 1
                    error('AnalyzeData: please only provide the channel number to be analyzed')
                end
                obj.channel = varargin{1};
            elseif strcmp( routine, 'all')
                obj.routine = routine;
            else
                error('AnalyzeData: unknown routine')
            end

            obj.path = path;

        end
        % }}}

        function loadFrame( path_frame)
            % given a path with .mat files of a single path, it loads them

            % look for files ending with 'global.mat', 'local.mat', 'globum.mat'

            % Get all file names that end in .mat
            path_frame = what( path_frame); path_frame = pathframe.path;
            files = dir( [path_frame, filesep, '*.mat'] );

        end

    end
