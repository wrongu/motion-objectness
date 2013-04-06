structGT_MOT = struct('annotations', {}, 'annot_dim', [], 'vid', '');

for c=1:length(classes)
    for i=1:10
        ind = (c-1)*10 + i;
        a = load(sprintf('Data/yourData/Youtube/annot/%s%.2d.mat', classes(c).name, i));
        f = floor(rand * a.num_frames)+1;
        structGT_MOT(ind).frame = f;
        annot = a.annotations{f};
        structGT_MOT(ind).boxes = ...
            [annot.xtl annot.ytl annot.xbr annot.ybr] + 1;
        structGT_MOT(ind).annot_dim = [a.width a.height];
        structGT_MOT(ind).vid = sprintf('Training/Videos/%s%.2d.mat', classes(c).name, i);
    end
end