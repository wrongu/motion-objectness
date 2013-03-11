% This function extracts the segmented point information from the Tracks
% structure at a given frame. See readTracksFile for definition of Tracks
% struct array; see ppmToSegments for definition of Segments struct

function S = sliceTracks(T, frame)
labels = [T.label];
n_segments = length(unique(labels));
S = struct('points', cell(1, n_segments));

for t=1:length(T)
    % find where, if at all, the trajectory T(t) goes through the given
    % frame
    i = T(t).frames == frame;
    S(T(t).label).points = horzcat(S(T(t).label).points, T(t).points(:,i));
end
end