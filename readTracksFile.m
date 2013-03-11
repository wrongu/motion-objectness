% This function reads the Tracks.dat file output by motion segmentation
%
% The return value Tracks is a struct array where each Track has fields
%   Tracks(i).label  % an integer label
%   Tracks(i).points % a 2xn matrix of points (as column vectors with (x,y))
%   Tracks(i).frames % a 1xn array containing the frame for each point
function Tracks = readTracksFile(filePath)
if ~exist(filePath, 'file')
    fprintf('Error: file %s does not exist\n', filePath);
elseif strfind(filePath, '.dat') ~= length(filePath)-3
    fprintf('Error: file must be .dat\n', filePath);
end

Tracks = struct('label', 0, 'points', []);

fid = fopen(filePath);
n_frames = fscanf(fid, '%d\n', 1); % first line is # frames
n_tracks = fscanf(fid, '%d\n', 1); % second line is # tracks total

% read each track (revers loop for quick allocation of Tracks)
for t=n_tracks:-1:1
    % new track parsing: first line is "<label> <n_pts>"
    Tracks(t).label = fscanf(fid, '%d ', 1);
    n_pts = fscanf(fid, '%d\n', 1);
    Tracks(t).points = zeros(2, n_pts);
    Tracks(t).frames = zeros(1, n_pts);
    % get each point
    for n=1:n_pts
        x = fscanf(fid, '%f ', 1);
        y = fscanf(fid, '%f ', 1);
        f = fscanf(fid, '%d\n', 1);
        Tracks(t).points(:,n) = [x; y];
        Tracks(t).frames(n) = f;
    end
end

fclose(fid);

% Change labels to be in 1:num_label_types
labels_ = [Tracks.label];
labels = zeros(size(labels_));
types = unique(labels_);
k=1;
for l=types
    labels(labels_==l) = k;
    k = k+1;
end
for t=1:n_tracks
    Tracks(t).label = labels(t);
end

end