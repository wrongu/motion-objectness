% read and annotate all videos

top_dir = fullfile('Data', 'yourData', 'Youtube');
classes = dir(fullfile(top_dir, '10class'));

for c=1:length(classes)
    class = classes(c);
    if class.isdir && class.name(1) ~= '.'
        vids = dir(fullfile(top_dir, '10class', class.name, '*.avi'));
        for v=1:length(vids)
            fprintf('%s%s\n', class.name, vids(v).name);
            vid = VideoReader(fullfile(top_dir, '10class', class.name, vids(v).name));
            a = load(fullfile(top_dir, 'annot', ...
                sprintf('%s%.2d.mat', class.name, v)));
            mov = box_video(vid, a.annotations);
            writer = VideoWriter(fullfile(top_dir, 'combined', ...
                sprintf('%s%.2d_annot.avi', class.name, v)));
        end
    end
end