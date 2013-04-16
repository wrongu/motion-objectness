function boxes = computeScoresWithMOT(descriptorGT,cue,params,windows)
V = VideoReader([params.trainingImages descriptorGT.vid]);
img = read(V, descriptorGT.frame);
vdata = get(V, {'NumberOfFrames'});
if nargin<4
    %no windows provided - so generate them -> single cues
    
    switch cue
        
        case 'MS' %Multi-scale Saliency
            
            xmin = [];
            ymin = [];
            xmax = [];
            ymax = [];
            score = [];
            img = gray2rgb(img); %always have 3 channels
            [height width ~] = size(img);
            
            for sid = 1:length(params.MS.scale) %looping over the scales
                
                scale = params.MS.scale(sid);
                threshold = params.MS.theta(sid);
                min_width = max(2,round(params.min_window_width * scale/width));
                min_height = max(2,round(params.min_window_height * scale/height));
                
                samples = round(params.distribution_windows/(length(params.MS.scale)*3)); %number of samples per channel to be generated
                
                for channel = 1:3 %looping over the channels
                    
                    saliencyMAP = saliencyMapChannel(img,channel,params.MS.filtersize,scale);%compute the saliency map - for the current scale & channel
                    
                    thrmap = saliencyMAP >= threshold;
                    salmap = saliencyMAP .* thrmap;
                    thrmapIntegralImage = computeIntegralImage(thrmap);
                    salmapIntegralImage =  computeIntegralImage(salmap);
                    
                    scoreScale = slidingWindowComputeScore(double(saliencyMAP), scale, min_width, min_height, threshold, salmapIntegralImage, thrmapIntegralImage);%compute all the windows score
                    %keyboard;
                    indexPositives = find(scoreScale>0); %find the index of the windows with positive(>0) score
                    scoreScale = scoreScale(indexPositives);
                    
                    indexSamples = scoreSampling(scoreScale, samples, 1);%sample from the distribution of the scores
                    scoreScale = scoreScale(indexSamples);
                    
                    [xminScale yminScale xmaxScale ymaxScale] = retrieveCoordinates(indexPositives(indexSamples) - 1,scale);%
                    xminScale = xminScale*width/scale;
                    xmaxScale = xmaxScale*width/scale;
                    yminScale = yminScale*height/scale;
                    ymaxScale = ymaxScale*height/scale;
                    
                    score = [score;scoreScale];
                    xmin = [xmin ;xminScale];
                    ymin = [ymin ;yminScale];
                    xmax = [xmax ;xmaxScale];
                    ymax = [ymax ;ymaxScale];
                    
                end%loop channel
                
            end%loop sid
            
            boxes = [xmin ymin xmax ymax score];
            boxes = boxes(1:params.distribution_windows,:);%might be more than 100.000
            
        case 'CC'
            
            windows = generateWindows(img, 'uniform', params);%generate windows
            boxes = computeScoresWithMOT(descriptorGT, cue, params, windows);
            
        case 'ED'
            
            windows = generateWindows(img, 'dense', params, cue);%generate windows
            boxes = computeScoresWithMOT(descriptorGT, cue, params, windows);
            
        case 'SS'
            windows = generateWindows(img,'dense', params, cue);
            boxes = computeScoresWithMOT(descriptorGT, cue, params, windows);
            
        case 'MOT'
            windows = generateWindows(img,'dense', params, cue);
            boxes = computeScoresWithMOT(descriptorGT, cue, params, windows);
    end
    
else
    %windows are provided so score them
    switch cue
        
        case 'CC'
            
            [height width ~] = size(img);
            
            imgLAB = rgb2lab(img);%get the img in LAB space
            Q = computeQuantMatrix(imgLAB,params.CC.quant);
            integralHistogram = computeIntegralHistogramMex(double(Q), height, width, prod(params.CC.quant));
            
            xmin = round(windows(:,1));
            ymin = round(windows(:,2));
            xmax = round(windows(:,3));
            ymax = round(windows(:,4));
            
            score = computeScoreContrast(double(integralHistogram), height, width, xmin, ymin, xmax, ymax, params.CC.theta, prod(params.CC.quant), size(windows,1));%compute the CC score for the windows
            boxes = [windows score];
            
        case 'ED'
            
            [~, ~, temp] = size(img);
            if temp==3
                edgeMap = edge(rgb2gray(img),'canny');%compute the canny map for 3 channel images
            else
                edgeMap = edge(img,'canny');%compute the canny map for grey images
            end
            
            h = computeIntegralImage(edgeMap);
            
            xmin = round(windows(:,1));
            ymin = round(windows(:,2));
            xmax = round(windows(:,3));
            ymax = round(windows(:,4));
            
            xmaxInner = round((xmax*(200+params.ED.theta)/(params.ED.theta+100) + xmin*params.ED.theta/(params.ED.theta+100)+100/(params.ED.theta+100)-1)/2);
            xminInner  = round(xmax + xmin - xmaxInner);
            ymaxInner = round((ymax*(200+params.ED.theta)/(params.ED.theta+100) + ymin*params.ED.theta/(params.ED.theta+100)+100/(params.ED.theta+100)-1) /2);
            yminInner  = round(ymax + ymin - ymaxInner);
            
            scoreWindows = computeIntegralImageScores(h,[xmin ymin xmax ymax]);
            scoreInnerWindows = computeIntegralImageScores(h,[xminInner yminInner xmaxInner ymaxInner]);
            areaWindows = (xmax - xmin + 1) .* (ymax - ymin +1);
            areaInnerWindows = (xmaxInner - xminInner + 1) .* (ymaxInner - yminInner + 1);
            areaDiff = areaWindows - areaInnerWindows;
            areaDiff(areaDiff == 0) = inf;
            
            score = ((xmax - xmaxInner + ymax - ymaxInner)/2) .* (scoreWindows - scoreInnerWindows) ./ areaDiff;
            boxes = [windows score];
            
        case 'SS'
            
            currentDir = pwd;
            soft_dir = params.SS.soft_dir;
            basis_sigma = params.SS.basis_sigma;
            basis_k = params.SS.theta;
            basis_min_area = params.SS.basis_min_area;
            imgType = params.imageType;
            imgBase = tempname(params.tempdir);%find a unique name for a file in params.tempdir
            imgBase = imgBase(length(params.tempdir)+1:end);
            %imgName = [imgBase '.' imgType];
            imgName = [imgBase '.ppm'];
            cd(params.tempdir);
            fprintf('imwrite(img@[%d x %d], "%s", "%s")\n', size(img, 2), size(img, 1), imgName, 'PPM');
            imwrite(img,imgName,'PPM');
            if ~exist(imgName, 'file')
                error('imwrite didnt work..?\n');
            end
            segmFileName = [imgBase '_segm.ppm'];
            
            if not(exist(segmFileName,'file'))
                % convert image to ppm
                %                 if not(strcmp(imgName(end-2:end), 'ppm'))
                %                     %I = imread(imgName);
                %                     imwrite(img, [imgBase '.ppm'], 'PPM');
                %                     %                     [~, convert_cmd] = system('which convert');
                %                     %                     cmd = [ convert_cmd ' "' imgName '" "' imgBase '.ppm"' ];
                %                     %                     system(cmd);
                %                     %                     clear convert_cmd;
                %                 end
                % setting segmentation params
                %                 I = imread([imgBase '.ppm']);
                Iarea = size(img,1)*size(img,2);
                sf = sqrt(Iarea/(300*200));
                sigma = basis_sigma*sf;
                min_area = basis_min_area*sf;
                k = basis_k;
                % segment image
                cmd = [soft_dir '/segment ' num2str(sigma) ' ' num2str(k) ' ' num2str(min_area) ' "' imgBase '.ppm' '" "' segmFileName '"' ];
                fprintf('running cmd:\n\t%s\n',cmd);
                system(cmd);
                % delete image ppm
                if not(strcmp(imgType, 'ppm'))
                    delete([imgBase '.ppm']);
                    delete(imgName);
                end
                % S is a 3-channel _image_ (loaded from ppm) with
                %   superpixels colored together
                S = imread(segmFileName);
                delete(segmFileName);
            else    % segmentation file found
                S = imread(segmFileName);
            end
            
            cd(currentDir);
            % convert 3-channel S to 1-channel N where N(i,j)=k means that
            % pixel (i,j) is in the kth cluster
            N = numerizeLabels(S);
            
            %             subplot(1,2,1);
            %             image(N*64/max(max(N)));
            %             set(gca, 'YDir', 'reverse');
            %             subplot(1,2,2);
            %             image(img);
            %             pause;
            
            % get full set of [c;r] pixel coords in each super pixel, and
            % area of each. IE superpixels(k).coords; and
            % superpixels(k).area
            superpixels = segmentArea(N);
            % break N apart and compute integral images (w+1, h+1, n_SS)
            %   dimensions
            integralHist = integralHistSuperpixels(N);
            
            xmin = round(windows(:,1));
            ymin = round(windows(:,2));
            xmax = round(windows(:,3));
            ymax = round(windows(:,4));
            
            areaSuperpixels = [superpixels(:).area];
            areaWindows = (xmax - xmin + 1) .* (ymax - ymin + 1);
            % size: num windows X num superpixels
            intersectionSuperpixels = zeros(length(xmin),size(integralHist,3));
            
            for dim = 1:size(integralHist,3)
                intersectionSuperpixels(:,dim) = computeIntegralImageScores(integralHist(:,:,dim),windows);
            end
            
            score = ones(size(windows,1),1) - (sum(min(intersectionSuperpixels,repmat(areaSuperpixels,size(windows,1),1) - intersectionSuperpixels),2)./areaWindows);
            boxes = [windows score];
            
            % New case for motion segmentation (same sort of superpixel
            %   straddling, but with segmented trajectory points)
            %   Most code copied from 'SS' case, above
        case 'MOT'
            % execute motion segmentation algorithm. Defaults will save
            % results to moseg2012/marple2/OchsBroxResults/
            % (max and min ensure the given frame is included, and their position
            % is chosen so that GT.frame is at the end of the sequence considered.
            sf = max(1, descriptorGT.frame - params.MOT.nframes + 1);
            ef = min(sf + params.MOT.nframes - 1, vdata{1});
            % the only parameter we can really control is sampling. Ideally
            % we would change the thresholds for segmentation, but we can't
            % access that aspect of the executable yet
            sampling = params.MOT.theta;
            nameparts = regexp(descriptorGT.vid, '[^.]+', 'match');
            class = nameparts{1};
            tracks_f = fullfile(params.bmf_locations, class, 'OchsBroxResults', ...
                ['Tracks' num2str(sf) '_' num2str(ef) '_' num2str(sampling) '.dat']);
            % recompute motion segmentation iff result file does not
            % exist
            if ~exist(tracks_f, 'file')
                % ensure the BMF file exists. If not, make it, use it,
                % delete it (assume it didnt exist already because of disk
                % space limitations)
                bmf_file = fullfile(params.bmf_locations, class, [class '.bmf']);
                delete_bmf = false;
                if ~exist(bmf_file, 'file')
                    fprintf('creating bmf version\n');
                    videoToBMF(fullfile(params.trainingImages, descriptorGT.vid), ...
                        descriptorGT.vid(1:end-4), fileparts(bmf_file));
                    delete_bmf = true;
                end
                ld_lib_cmd = ['export LD_LIBRARY_PATH=~/lib:' ...
                    '~/motion-objectness/moseg2012:$LD_LIBRARY_PATH'];
                cmd = [params.MOT.executable ' ' bmf_file ...
                    ' ' num2str(sf) ' ' num2str(ef)  ' ' ...
                    num2str(sampling)];
                fprintf('%s\n------------\n', cmd);
                tstart = tic;
                system([ld_lib_cmd '; ' cmd]);
                fprintf('motion segmentation algorithm finished in %d seconds\n', toc(tstart));
                outfile = fullfile(params.bmf_locations, class, 'OchsBroxResults', ...
                    ['Tracks' num2str(ef-sf+1) '.dat']);
                movefile(outfile, tracks_f);
                if delete_bmf
                    fprintf('deleting bmf file\n');
                    delete(bmf_file);
                end
            end
            % load trajectories from file
            fprintf('reading tracks file %s\n', tracks_f);
            Tracks = readTracksFile(tracks_f);
            % slice out the frame we want
            disp('slicing');
            S = sliceTracks(Tracks, descriptorGT.frame);
            
            % sanity-check plot:
            %                 colors = 'rbgcyk';
            %                 h = figure();
            %                 hold on;
            %                 for s=1:length(S)
            %                     scatter(S(s).points(1,:), S(s).points(2,:), colors(s));
            %                 end
            %                 hold off;
            %                 set(gca, 'YDir', 'reverse');
            %                 pause;
            %                 close(h);
            
            N = sliceToSuperPixels(S, size(img, 2), size(img,1));
            % The following is the same as the above  (case 'SS')
            
            % superpixels is a struct array with 'points' and 'area'
            % fields
            superpixels = segmentArea(N);
            % integralHist is a 3D matrix where the kth layer is the
            % integral image for the kth superpixel.
            %   size(integralHist,3) := num of superpixels
            integralHist = integralHistSuperpixels(N);
            
            xmin = round(windows(:,1));
            ymin = round(windows(:,2));
            xmax = round(windows(:,3));
            ymax = round(windows(:,4));
            
            areaSuperpixels = [superpixels(:).area];
            areaWindows = (xmax - xmin + 1) .* (ymax - ymin + 1);
            
            intersectionSuperpixels = zeros(length(xmin),size(integralHist,3));
            
            for dim = 1:size(integralHist,3)
                % compute the sum of the integral image _within each
                % window_. Since the IIs are binary, this is the
                % _count_ of the number of pixels from each SP in each
                % window.
                % Each ROW corresponds to a window. each COL to a
                % superpixel.
                intersectionSuperpixels(:,dim) = computeIntegralImageScores(integralHist(:,:,dim),windows);
            end
            % according to the paper:
            % score = 1 - [SUM over superpixels](min([area inside window], [are outside window]) / [area of window])
            score = ones(size(windows,1),1) - ...
                (sum( ...
                min( ... % note that min(), here, is element-wise since the two arguments are matrices of the same size
                intersectionSuperpixels, ...
                repmat(areaSuperpixels,size(windows,1), 1) - intersectionSuperpixels ...
                ), ...
                2)./areaWindows);
            boxes = [windows score];
            
            
        otherwise
            error('Option not known: check the cue names');
    end
end

end



function saliencyMAP = saliencyMapChannel(inImg,channel,filtersize,scale)

inImg = im2double(inImg(:,:,channel));
inImg = imresize(inImg,[scale,scale],'bilinear');

%Spectral Residual
myFFT = fft2(inImg);
myLogAmplitude = log(abs(myFFT));
myPhase = angle(myFFT);
mySmooth = imfilter(myLogAmplitude,fspecial('average',filtersize),'replicate');
mySpectralResidual = myLogAmplitude-mySmooth;
saliencyMAP = abs(ifft2(exp(mySpectralResidual+1i*myPhase))).^2;

%After Effect
saliencyMAP = imfilter(saliencyMAP,fspecial('disk',filtersize));
saliencyMAP = mat2gray(saliencyMAP);
end
