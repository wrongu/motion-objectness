% read and annotate all videos

top_dir = fullfile('Data', 'yourData', 'Youtube');
classes = dir(fullfile(top_dir, '10class'));

for c=1:length(classes)
    class = classes(c);
    cname = class.name;
    if class.isdir && cname(1) ~= '.'
        vids = dir(fullfile(top_dir, '10class', cname, '*.avi'));
        parfor v=1:length(vids)
            fprintf('%s%s\n', cname, vids(v).name);
            if ~exist(fullfile(top_dir, 'combined', ...
                    sprintf('%s%.2d_annot.avi', cname, v)), 'file')
                vid = VideoReader(fullfile(top_dir, '10class', cname, vids(v).name));
                a = load(fullfile(top_dir, 'annot', ...
                    sprintf('%s%.2d.mat', cname, v)));
                mov = box_video(vid, a.annotations);
                writer = VideoWriter(fullfile(top_dir, 'combined', ...
                    sprintf('%s%.2d_annot.avi', cname, v)));
                open(writer);
		writeVideo(writer, mov);
            end
        end
    end
end
