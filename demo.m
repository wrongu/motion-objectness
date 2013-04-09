%% DEMO ORIGINAL
% imgExample = imread('002053.jpg');
% boxes = runObjectness(imgExample,10);
% figure,imshow(imgExample),drawBoxes(boxes);

%% DEMO MOTION

% params.MOT.bmfFile = 'moseg2012/marple2/marple2.bmf';

img.vid = 'bear06.avi';
img.frame = 12;
boxes = runMotionObjectness(img, 10,params);
V = VideoReader(fullfile(params.trainingImages, img.vid));
I = read(V, img.frame);
figure;
imshow(I);
drawBoxes(boxes);