% This function reads the ppm results from Brox-Malik motion based
% segmentation and saves a .mat file with the movie results.
%
% Author: Richard Lange
% Date Created: March 10, 2013

function movies = movieResults(images_path, start, finish)
    if(~exist('readFlowFile', 'file'))
        if(exist('flow-code-matlab', 'dir'))
            addpath(fullfile(pwd, 'flow-code-matlab'));
        else
            fprintf('movieResults: must have the directory flow-code-matlab in the search path\n');
        end
        return;
    end
    % prepare platform-independent file name format strings
    ff = safe_fullfile('images_path', 'BroxMalikResults', 'ForwardFlow%s.flo');
    bf = safe_fullfile('images_path', 'BroxMalikResults', 'BackwardFlow%s.flo'); 
    tr = safe_fullfile('images_path', 'BroxMalikResults', 'Tracking%s.ppm'); 
    sm = safe_fullfile('images_path', 'BroxMalikResults', 'Segments%s.ppm');
    im = 'marple2_%s.ppm';
    
    movies(1).file = im; movies(1).movie = empty_movie();
    movies(2).file = ff; movies(2).movie = empty_movie();
    movies(3).file = bf; movies(3).movie = empty_movie();
    movies(4).file = tr; movies(4).movie = empty_movie();
    movies(5).file = sm; movies(5).movie = empty_movie();
    
    for m=1:5
        disp(movies(m).file);
        frame_i = 1;
        for frame=start:finish
            s = num2str(frame);
            while length(s) < 3
                s = ['0' s]; % because sprintf('%3d', frame) isn't working
            end
            final_path = safe_fullfile(images_path, sprintf(movies(m).file, s));
            if(~exist(final_path, 'file'))
                fprintf('%s does not exist');
                continue;
            end
            disp(final_path);
            if strfind(movies(m).file, 'flo')
                I = readFlowFile(final_path);
                I = flowToColor(I);
            else
                I = imread(final_path);
            end
            imshow(I);
            movies(m).movie(frame_i) = getframe();
            frame_i = frame_i + 1;
        end
    end
    
    % save results
    save(fullfile(images_path, 'movieResults.mat'), 'movies');
end

function m = empty_movie()
    m = struct('cdata', [], 'colormap', []);
end

function s = safe_fullfile(varargin)
    % replace windows '\' with escaped '\\'
    % the regular expression prevents further modifications like
    % '\\' becoming '\\\\'
    s = regexprep(fullfile(varargin{:}), '(\w)\\(\w)', '$1\\\\$2');
end