classdef ImageProcessingProcess < Process
    %A class definition for a generic image processing process. That is, a
    %process which takes in images and produces images of the same
    %dimension and number as output. These images may or may not overwrite
    %the original input images.
    %
    %
    % Hunter Elliott, 5/2010
    %
    
    methods (Access = public)
        
        function obj = ImageProcessingProcess(owner,name,funName,funParams,...
                                              inImagePaths,outImagePaths)
                                          
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;                
            end
            
            obj = obj@Process(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;                              
            end
            if nargin > 3
               obj.funParams_ = funParams;              
            end
            
            if nargin > 4
              if ~isempty(inImagePaths) && numel(inImagePaths) ...
                      ~= numel(owner.channels_) || ~iscell(inImagePaths)
                 error('lccb:set:fatal','Input image paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end         
                obj.inFilePaths_ = inImagePaths;
            else
                %Default is to use raw images as input.
                obj.inFilePaths_ = owner.getChannelPaths;               
            end                        
            if nargin > 5               
                if ~isempty(outImagePaths) && numel(outImagePaths) ... 
                        ~= numel(owner.channels_) || ~iscell(outImagePaths)
                    error('lccb:set:fatal','Output image paths must be a cell-array of the same size as the number of image channels!\n\n'); 
                end
                obj.outFilePaths_ = outImagePaths;              
            else
                obj.outFilePaths_ = cell(1,numel(owner.channels_));               
            end
            
        end
        
        function setOutImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n'); 
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
               if ~exist(imagePath{j},'dir')
                   error('lccb:set:fatal',...
                       ['The directory specified for channel ' ...
                       num2str(chanNum(j)) ' is invalid!']) 
               else
                   obj.outFilePaths_{1,chanNum(j)} = imagePath{j};                
               end
            end
            
            
        end
        function clearOutImagePath(obj,chanNum)

            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n'); 
            end
            
            for j = 1:numel(chanNum)
                obj.outFilePaths_{1,chanNum(j)} = [];
            end
        end
        function setInImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n'); 
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            
            for j = 1:nChan
                if ~obj.owner_.isBF
                   if ~exist(imagePath{j},'dir')
                       error('lccb:set:fatal',...
                           ['The directory specified for channel ' ...
                           num2str(chanNum(j)) ' is invalid!']) 
                   else
                       obj.inFilePaths_{1,chanNum(j)} = imagePath{j};                
                   end                
                else
                    if ~exist(imagePath{j},'file')
                       error('lccb:set:fatal',...
                           ['The file specified for channel ' ...
                           num2str(chanNum(j)) ' is invalid!']) 
                    else
                       obj.inFilePaths_{1,chanNum(j)} = imagePath{j};                
                    end     
                end
            end                        
        end
        function fileNames = getOutImageFileNames(obj,iChan)
            nChanTot = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan), iChan = 1:nChanTot; end
            if obj.checkChannelOutput(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.outFilePaths_(1,iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)                    
                    error('Incorrect number of images found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end    
            
            
        end
        function fileNames = getInImageFileNames(obj,iChan)
            nChanTot = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan), iChan = 1:nChanTot; end
           
            if obj.checkChanNum(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.inFilePaths_(1,iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)                    
                    error('Incorrect number of images found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end    
            
            
        end
        
        function status = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channels have valid output images          
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan), iChan = 1:nChanTot; end
           assert(all(obj.checkChanNum(iChan)));
           status =  arrayfun(@(x) exist(obj.outFilePaths_{1,x},'dir') && ...
               ~isempty(imDir(obj.outFilePaths_{1,x})),iChan);
        end
        
        
        
        function outIm = loadOutImage(obj,iChan,iFrame)
            outIm=obj.loadChannelOutput(iChan,iFrame);
        end
        
        function outIm = loadChannelOutput(obj,iChan,iFrame,varargin)
             % Input check
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @obj.checkChanNum);
            ip.addRequired('iFrame', @obj.checkFrameNum);
            if obj.owner_.is3D()
                ip.addOptional('iZ', @obj.checkDepthNum);
            end
            ip.addParamValue('output',[],@ischar);            
            ip.parse(obj,iChan,iFrame,varargin{:})
            imNames = obj.getOutImageFileNames(iChan);
            if obj.getOwner().is3D()
                %iZ = ip.Results.iZ;
                if ~exist('iZ', 'var')
                    outIm = tif3Dread([obj.outFilePaths_{1,iChan} filesep imNames{1}{iFrame}]);
                else 
                    outIm =imread([obj.outFilePaths_{1,iChan} filesep imNames{1}{iFrame}], iZ);
                end
            else
                outIm =imread([obj.outFilePaths_{1,iChan} filesep imNames{1}{iFrame}]);
            end
        end
        
        function drawImaris(obj,iceConn)
            
            dataSet = iceConn.mImarisApplication.GetDataSet;            
            nChanRaw = numel(obj.owner_.channels_);
            nFrames = obj.owner_.nFrames_;
            for iChan = 1:nChanRaw
                
                if obj.checkChannelOutput(iChan)
                                        
                    nChanDisp = dataSet.GetSizeC + 1;
                    dataSet.SetSizeC(nChanDisp)
                    dataSet.SetChannelName(nChanDisp-1,[ char(dataSet.GetChannelName(iChan-1)) ' ' obj.name_]);
                    dataSet.SetChannelColorRGBA(nChanDisp-1,dataSet.GetChannelColorRGBA(iChan-1));                    
                    datMin = Inf;
                    datMax = -Inf;
                    for iFrame = 1:nFrames
                        vol = obj.loadChannelOutput(iChan,iFrame);                        
                        datMin = min(min(vol(:)),datMin);
                        datMax = max(max(vol(:)),datMax);
                        iceConn.setDataVolume(vol,nChanDisp-1,iFrame-1);                       
                    end
                    dataSet.SetChannelRange(nChanDisp-1,datMin,datMax);                   
                end                
            end
            
        end
        
        
    end
    methods(Static)
        function output = getDrawableOutput()
            output(1).name='Images';
            output(1).var='';
            output(1).formatData=@mat2gray;
            output(1).type='image';
            output(1).defaultDisplayMethod=@ImageDisplay;
        end
        
    end
    
end