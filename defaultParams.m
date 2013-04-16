function params = defaultParams(dirRoot)

if nargin < 1
    dirRoot = pwd;
end

%params in general
params.min_window_height = 10;
params.min_window_width  = 10;
params.distribution_windows = 100000;
params.sampled_windows = 1000;
params.trainingImages = [dirRoot '/Training/Videos/'];
params.trainingLabels = [dirRoot '/Data/yourData/Youtube/annot/'];
params.trainingExamples = [dirRoot '/Training/Videos/Examples/'];
params.bmf_locations = [dirRoot '/Training/BMF/'];
% params.trainingImages = [dirRoot '/Training/Images/'];
% params.trainingExamples = [dirRoot '/Training/Images/Examples/'];
params.imageType = 'JPG';
params.data = [dirRoot '/Data/'];
params.yourData = [dirRoot '/Data/yourData/'];
params.pobj = 0.0797;
params.tempdir = [dirRoot '/tmpdir/'];
params.pascalThreshold = 0.5;
params.cues = {'MS','CC','SS','MOT'};%full objectness measure
params.sampling = 'nms';%alternative sampling method - 'multinomial'

%params for MS
params.MS.name = 'Multiscale-Saliency';
params.MS.colortype = 'rgb';
params.MS.filtersize = 3;
params.MS.scale = [16 24 32 48 64];
params.MS.theta = [0.43 0.32 0.34 0.35 0.26];
params.MS.domain = repmat(0.01:0.01:1,5,1);
params.MS.sizeNeighborhood = 7;
params.MS.bincenters = 0:1:500;
params.MS.numberBins = length(params.MS.bincenters) - 1;

%params for CC
params.CC.name = 'Color-Contrast';
params.CC.theta = 100;
params.CC.domain = 100:1:200;
params.CC.quant = [4 8 8];
params.CC.bincenters = 0:0.01:2;
params.CC.numberBins = length(params.CC.bincenters) - 1;

%params for ED
params.ED.name = 'Edge-Density';
params.ED.theta = 17;
params.ED.domain = 1:1:100;
params.ED.crop_size = 200;
params.ED.pixelDistance = 8;
params.ED.imageBorder = 0;
params.ED.bincenters = 0:0.05:5;
params.ED.numberBins = length(params.ED.bincenters) - 1;

%params for SS
params.SS.name = 'Superpixels-Straddling';
params.SS.basis_sigma = 0.5;
params.SS.theta = 450;
params.SS.domain = 200:25:2000;
params.SS.basis_min_area = 200;
params.SS.soft_dir =[dirRoot '/segment/'];
params.SS.pixelDistance = 8;
params.SS.imageBorder = 0.05;
params.SS.bincenters = 0:0.01:1;
params.SS.numberBins = length(params.SS.bincenters) - 1;

%parmas for MOT (Motion segmentation)
params.MOT.name = 'Motion Segmentation';
params.MOT.theta = 8;
params.MOT.pixelDistance = 8;
params.MOT.imageBorder = 0.05
params.MOT.nframes = 50;
params.MOT.sampling = 8;
params.MOT.domain = [1 2 4 6 8 10 12 14 16 18 20];
params.MOT.bincenters = 0:0.01:1;
params.MOT.executable = fullfile(dirRoot, 'moseg2012', 'motionsegOB');
% params.MOT.bmfFile = fullfile(dirRoot, 'moseg', 'marple2', 'marple2.bmf');
% params.MOT.executable = fullfile(dirRoot, 'moseg', 'motionsegBM');
% params.MOT.resultsDir = fullfile(dirRoot, 'moseg', 'marple2', 'BroxMalikResults');
% params.MOT.bmfFile = fullfile(dirRoot, 'moseg2012', 'marple2', 'marple2.bmf');
% params.MOT.resultsDir = fullfile(dirRoot, 'moseg2012', 'marple2', 'OchsBroxResults');
