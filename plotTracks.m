% This function plots the struct array of Tracks (see readTracksFile.m)
%   in 3D where x,y are the horizontal and veritcal image axes, and
%   z is the time axis (frames)

function h = plotTracks(T, h_in)
if nargin < 2
    h = figure();
else
    h = h_in;
end

n_frames = 0;
n_tracks = length(T);
for i=1:n_tracks
    n_frames = max(n_frames, length(T(i).frames));
end

T_mat = zeros(3, n_frames, n_tracks);
for i=1:n_tracks
    T_mat(:, 1:length(T(i).frames), i) = [T(i).points; T(i).frames];
end

labels = [T.label];
types = unique(labels);

colors = hsv(length(types));

figure(h);
hold on;
for i=1:n_tracks
    range = 1:length(T(i).frames);
    % y/z swapped for better visibility
    plot3(T_mat(1, range, i), T_mat(3, range, i), T_mat(2, range, i), ...
        'Color', colors(labels(i),:));
end
hold off;
xlabel('x');
ylabel('frame');
zlabel('y');
set(gca, 'YTick', 1:n_frames);
set(gca, 'YDir', 'reverse');
end