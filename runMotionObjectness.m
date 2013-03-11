% Run the combined objectness and motion segmentation algorithm from
% scratch (including learning parameters)
%
% Author: Richard Lange
% Date Created: March 11, 2013

function boxes = runMotionObjectness(numberSamples, params)
dir_root = pwd;%change this to an absolute path
img = params.MOT.frame;
if nargin < 2
    try            
        struct = load([dir_root '/Data/params.mat']);
        params = struct.params;
        clear struct;
    catch
        params = defaultParams(dir_root);
        save([dir_root '/Data/params.mat'],'params');
    end
end

if length(params.cues)==1    
    %single cues
    
    distributionBoxes = computeScoresWithMOT(img,params.cues{1},params); 
    
    switch lower(params.sampling)
        case 'nms'
            %nms sampling
            
            %consider only params.distribution_windows (= 100k windows)
            if size(distributionBoxes,1) > params.distribution_windows
                indexSamples = scoreSampling(distributionBoxes(:,5),params.distribution_windows,1);
                distributionBoxes = distributionBoxes(indexSamples,:);
            end
            
            %sampling
            boxes = nms_pascal(distributionBoxes, 0.5,numberSamples);
            
        case 'multinomial'            
            %multinomial sampling
            
            %sample from the distribution of the scores
            indexSamples = scoreSampling(distributionBoxes(:,end),numberSamples,1);
            boxes = distributionBoxes(indexSamples,:);
                        
        otherwise
            display('sampling procedure unknown')
    end

else
    %combination of cues
    
    if not(ismember('MS',params.cues)) 
        display('ERROR: combinations have to include MS');
        boxes = [];
        return
    end
    
    if length(unique(params.cues)) ~= length(params.cues)
        display('ERROR: repetead cues in the combination');
        boxes = [];
        return
    end
    
    distributionBoxes = computeScoresWithMOTWithMOT(img,'MS',params);    
    %rearrange the cues such that 'MS' is the first cue
    if ~strcmp(params.cues{1},'MS')
        params.cues{strcmp(params.cues,'MS')} = params.cues{1};
        params.cues{1} ='MS';
    end
    
    score = zeros(size(distributionBoxes,1),length(params.cues));
    score(:,1) = distributionBoxes(:,end);
    windows = distributionBoxes(:,1:4);
    for idx = 2:length(params.cues)
        temp = computeScoresWithMOT(img,params.cues{idx},params,windows);
        score(:,idx) = temp(:,end);       
    end
    scoreBayes = integrateBayes(params.cues,score,params);      
    
    switch lower(params.sampling)
        case 'nms'
            %nms sampling
                        
            distributionBoxes(:,5) = scoreBayes;
            boxes = nms_pascal(distributionBoxes, 0.5, numberSamples);
            
        case 'multinomial'            
            %multinomial sampling
            
            %sample from the distribution of the scores
            indexSamples = scoreSampling(scoreBayes,numberSamples,1);  
        boxes = [windows(indexSamples,:) scoreBayes(indexSamples,:)];   
                        
        otherwise
            display('sampling procedure unknown')
    end  
end