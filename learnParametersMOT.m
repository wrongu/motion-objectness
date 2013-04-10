% TODO - redefine the structGT thing so it works with annotated video data
% for training examples

function params = learnParametersMOT(pathNewTrainingFolder,dir_root,skip_MS)
%learns the parameters of the objectness function: theta_MS (for 5 scales),
%theta_CC, theta_ED, theta_SS, theta_MOT and also the likelihoods corresp to each cue

%dir_root - path where the software is installed - see README Setting things up
if nargin < 2
    dir_root = [pwd '/'];
end
params = defaultParams(dir_root);


if nargin < 3
    skip_MS = false;
end

if nargin == 1
    %train the parameters from another dataset
    params.trainingImages = pathNewTrainingFolder;
    origDir = pwd;
    cd(params.trainingImages);
    mkdir('Examples');
    params.trainingExamples = [pathNewTrainingFolder '/Examples/'];
    cd(origDir);
end

if ~skip_MS || ~exist([params.trainingExamples '/posnegMS.mat'], 'file')
    if skip_MS && ~exist([params.trainingExamples '/posnegMS.mat'], 'file')
        fprintf('computing MS because posnegMS.mat does not exist');
    end
    if skip_MS && exist([params.trainingExamples '/saveMSTheta.mat'], 'file')
        load([params.trainingExamples '/saveMSTheta.mat']);
        fprintf('skipped thetaMS calculation\n');
    else
        %learn parameters for MS
        for idx = 1: length(params.MS.scale)
            scale = params.MS.scale(idx);   
            params.MS.theta(idx) = learnThetaMSMOT(params,scale); 
        end

        save([params.trainingExamples '/saveMSTheta.mat'], 'params');
    end

    try 
        struct = load([params.trainingExamples '/posnegMS.mat'] );
        posnegMS = struct.posnegMS;
        clear struct;
        fprintf('loaded posnegMS\n');
    catch    
        fprintf('creating posnegMS\n');
        posnegMS = generatePosNegMS_MOT(params);
        save([params.trainingExamples '/posnegMS.mat'],'posnegMS');
    end
else
    load([params.trainingExamples '/posnegMS.mat']);
end

[likelihood, pObj] = deriveLikelihoodMS(posnegMS,params);
save([params.yourData 'MSlikelihood.mat'],'likelihood');
params.pObj = pObj;

%learn parameters for CC, ED, SS, MOT
cues = {'CC','ED','SS','MOT'};

for cid = 1:length(cues)
    cue = cues{cid};
    [thetaOpt, ~, ~] = learnThetaWithMOT(cue,params);
    params.(cue).theta = thetaOpt;    
    save([params.yourData upper(cue) 'likelihood'],'likelihood');
end
    
save([params.yourData '/params.mat'],'params');

end


function posneg = generatePosNegMS_MOT(params)

struct = load([params.trainingImages 'structGT_MOT.mat']);
structGT = struct.structGT_MOT;

for idx = length(structGT):-1:1
    V = VideoReader(fullfile(params.trainingImages, structGT(idx).vid));
    img = read(V, structGT(idx).frame);
    boxes = computeScoresWithMOT(structGT(idx),'MS',params);        
    posneg(idx).examples = boxes(:,1:4);
    labels = -ones(size(boxes,1),1);
    for idx_window = 1:size(boxes,1)        
        labelled_boxes = structGT(idx).boxes;
        % if the box from computeScores has a high enough pascal score with
        %   _any_ of the labels, score it 1 and break
        for bb_idx = 1:size(labelled_boxes,1)
            box = labelled_boxes(bb_idx,:);
            pascalScore = computePascalScore(box, boxes(idx_window,1:4));
            if (pascalScore >= params.pascalThreshold)
                labels(idx_window) = 1;
                break;
            end
        end
    end
    posneg(idx).labels = labels;
    posneg(idx).img = img;
    posneg(idx).scores = boxes(:,5);
end

end

function [likelihood pobj] = deriveLikelihoodMS(posneg,params)

examplesPos = zeros(length(posneg) * params.distribution_windows,1);
examplesNeg = zeros(length(posneg) * params.distribution_windows,1);

pos = 0; 
neg = 0;

for idx = 1:length(posneg)        
           
    indexPositive = find(posneg(idx).labels == 1);
    examplesPos(pos+1:pos+length(indexPositive)) = posneg(idx).scores(indexPositive);
    pos = pos + length(indexPositive);
    
    indexNegative = find(posneg(idx).labels == -1);
    examplesNeg(neg+1:neg+length(indexNegative)) = posneg(idx).scores(indexNegative);
    neg = neg + length(indexNegative);
           
end

examplesPos(pos+1:end) = [];
examplesNeg(neg+1:end) = [];

pobj = pos/(pos+neg);

posLikelihood = hist(examplesPos,params.MS.bincenters)/length(examplesPos) + eps;
negLikelihood = hist(examplesNeg,params.MS.bincenters)/length(examplesNeg) + eps;

likelihood = [posLikelihood;negLikelihood];
end
