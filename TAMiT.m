function TAMiT( varargin)
    % --------------------------------------------------------
    % --------------------------------------------------------
    % --------------------------------------------------------
    %     _______       __  __ _ _______ 
    %    |__   __|/\   |  \/  (_)__   __|
    %       | |  /  \  | \  / |_   | |   
    %       | | / /\ \ | |\/| | |  | |   
    %       | |/ ____ \| |  | | |  | |   
    %       |_/_/    \_\_|  |_|_|  |_|   
    % 
    % Toolkit for Automated Microtubule Tracking
    %
    % BSD 3-Clause License
    % 
    % Copyright (c) 2022 Saad J. Ansari
    %
    %   Optional Parameters:
    %       1. 'Display'    : Flag (0 or 1)
    %       2. 'LOC'        : String ('Local' or 'Summit')
    %
    %   Usage:
    %       TAMiT()             : Display OFF
    %       TAMiT('Display',1)  : Display ON
    %       
    %   Advanced Usage:
    %       TAMiT('LOC', 'Local')   : Local machine, 
    %       TAMiT('LOC', 'Summit')  : Summit supercomputer
    %
    % --------------------------------------------------------
    % --------------------------------------------------------
    % --------------------------------------------------------
    %
    % BSD 3-Clause License
    % 
    % Copyright (c) 2022 Saad J. Ansari
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are met:
    % 
    % 1. Redistributions of source code must retain the above copyright notice, this
    %    list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright notice,
    %    this list of conditions and the following disclaimer in the documentation
    %    and/or other materials provided with the distribution.
    % 
    % 3. Neither the name of the copyright holder nor the names of its
    %    contributors may be used to endorse or promote products derived from
    %    this software without specific prior written permission.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    % --------------------------------------------------------
    % --------------------------------------------------------
    
    % ----------------------------- PREP ---------------------------

    opts = parseArgs( varargin{:} );

    clc; close all;
    clearvars -except CFG opts

    % Make sure we have the correct paths
%     warning('off', 'MATLAB:rmpath:DirNotFound');
%     rmpath( genpath(pwd) );
%     warning('on', 'MATLAB:rmpath:DirNotFound');
    addpath( pwd);

    % check if params file exists in the current folder
    if exist( fullfile( pwd, 'initParams.m' ) ) ~= 2
        error( ['copy initParams.m to the current folder location : ', pwd] );
    end

    % run the settings file : creates a params.mat file in the save directory folder
    paramsPath = feval('initParams', opts);

    % Define cleanup tasks
    c2 = onCleanup( @() delete( gcp('nocreate') ) );
    c3 = onCleanup( @() disp('Closing files and cleaning up') );

    % ----------------------------- MAIN ---------------------------
   
    for jCell = 1 : length( paramsPath)

        % Run a single TAMiT Cell
        TAMiT_cell( paramsPath{ jCell} );

    end

    % ---------------------------- CLEANUP -------------------------
    
    % parseArgs {{{
    function opts = parseArgs( varargin)

        % default
        defaultDisplay = 0;
        defaultLOC = 'Local';

        % valid Arguments
        validDisplay = @(x) x==0 || x==1;
        validLOC = @(x) strcmpi( x, 'Local') || strcmpi( x, 'Summit') || strcmpi( x, 'Rumor');

        % Input Parser
        p = inputParser;
        addParameter( p, 'Display', defaultDisplay);
        addParameter( p, 'LOC', defaultLOC);

        parse( p, varargin{:});
        opts = p.Results;

    end
    % }}}
    
end
