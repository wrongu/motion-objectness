% NOTE - the end result of all these videos being computed adds up to about
% 5GB of space. For this reason, the code was changed to convert avi to bmf
% ad hoc, deleting it immediately after
vids = dir('Training/Videos/*.avi');
parfor i=1:100
    if ~exist(['Training/BMF/' vids(i).name(1:end-4)], 'dir') || ...
          isempty(dir(['Training/BMF/' vids(i).name(1:end-4) '/*.ppm']));
        fprintf('converting %s to bmf\n',vids(i).name);
        videoToBMF(['Training/Videos/' vids(i).name], vids(i).name(1:end-4), 'Training/BMF');
    end
end