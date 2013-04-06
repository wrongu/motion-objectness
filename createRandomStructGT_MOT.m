% this function creates a new set of random training examples (video)

function structGT_MOT = createRandomStructGT_MOT(n_ex, params)
    vids = dir([params.trainingImages '/*.avi']);
    vids = {vids.name};
    for idx=n_ex:-1:1
        % choose random video
        r = floor(rand*length(vids))+1;
        structGT_MOT(idx).vid = vids{r};
        % load annotation (.mat) file for this video
        datname = [structGT_MOT(idx).vid(1:end-3) 'mat'];
        a = load([params.trainingLabels datname]);
        % choose random frame
        f = floor(rand*a.num_frames)+1;
        structGT_MOT(idx).frame = f;
        % save box annotations for this frame
        annot = a.annotations{f};
        structGT_MOT(idx).boxes = ...
            [annot.xtl annot.ytl annot.xbr annot.ybr] + 1;
        structGT_MOT(idx).annot_dim = [a.width a.height];
    end
end