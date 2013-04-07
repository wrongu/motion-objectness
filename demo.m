%% DEMO ORIGINAL
% imgExample = imread('002053.jpg');
% boxes = runObjectness(imgExample,10);
% figure,imshow(imgExample),drawBoxes(boxes);

%% DEMO MOTION

% params.MOT.bmfFile = 'moseg2012/marple2/marple2.bmf';

boxes = runMotionObjectness(10,params);
figure,imshow(fullfile(fileparts(params.MOT.resultsDir), params.MOT.frame)),drawBoxes(boxes);