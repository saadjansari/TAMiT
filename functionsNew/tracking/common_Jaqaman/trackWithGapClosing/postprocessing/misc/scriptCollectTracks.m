
%% define batch job locations

%track locations
trackDir = {...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/Unstimulated/Field1/analysis/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/Unstimulated/Field2/analysis/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/Unstimulated/Field3/analysis/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/0min(A very minor stage drift occured while adding TSP1)/analysisBefore/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/0min(A very minor stage drift occured while adding TSP1)/analysisAfter/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/7min/analysis/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/14min/analysis/UTrackPackage/Tracking/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/28min/analysis/UTrackPackage/Tracking/',...
   };

%track file names
trackFile = {...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    'Channel_1_tracking_result.mat',...
    };

%diffusion analysis location
diffAnDir = {...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/Unstimulated/Field1/analysis/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/Unstimulated/Field2/analysis/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/Unstimulated/Field3/analysis/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/0min(A very minor stage drift occured while adding TSP1)/analysisBefore/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/0min(A very minor stage drift occured while adding TSP1)/analysisAfter/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/7min/analysis/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/14min/analysis/UTrackPackage/MotionAnalysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/DATA_SMI_LiveCellMovies_Revised_InformationCorrect/20120513_CD36_Lifeact/StimulatedTSP1/28min/analysis/UTrackPackage/MotionAnalysis/',...
    };

%diffusion file name
diffAnFile = {...
    'channel_1.mat',...
    'channel_1.mat',...
    'channel_1.mat',...
    'channel_1.mat',...
    'channel_1.mat'...
    'channel_1.mat',...
    'channel_1.mat'...
    'channel_1.mat',...
    };

%movie names
movieName = {...
    'unstim1',...
    'unstim2',...
    'unstim3',...
    '0min_before'...
    '0min_after'...
    '7min',...
    '14min',...
    '28min',...
    };

%calculate number of movies
numMovies1 = length(trackDir);
numMovies2 = length(trackFile);
numMovies3 = length(movieName);
if any([numMovies1-numMovies2 numMovies2-numMovies3 numMovies3-numMovies1]) ~= 0
    disp('Number of entries in movieInfoDir, movieInfoFile and saveResDir does not match!')
    return
end

%% collect movies

%generate structure for putting all tracks together
trackData = repmat(struct('name',[],'tracks',[],'diffAnalysis',[]),numMovies1,1);

for iMovie = 1 : numMovies1
    
    %read tracks
    load(fullfile(trackDir{iMovie},trackFile{iMovie}));
    load(fullfile(diffAnDir{iMovie},diffAnFile{iMovie}));
    
    %collect information
    trackData(iMovie).name = movieName{iMovie};
    trackData(iMovie).tracks = tracksFinal;
    trackData(iMovie).diffAnalysis = diffAnalysisRes;
    
end
