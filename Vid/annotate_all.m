% read and annotate all videos

if ~exist('force_recalc', 'var')
    force_recalc = false;
end

top_dir = fullfile('Data', 'yourData', 'Youtube');
classes = dir(fullfile(top_dir, '10class'));

for c=1:length(classes)
    class = classes(c);
    cname = class.name;
    if class.isdir && cname(1) ~= '.'
        vids = dir(fullfile(top_dir, '10class', cname, '*.avi'));
        parfor v=1:length(vids)
            fprintf('%s%s\n', cname, vids(v).name);
            if force_recalc || ~exist(fullfile(top_dir, 'combined', ...
                    sprintf('%s%.2d_annot.avi', cname, v)), 'file')
                vid = VideoReader(fullfile(top_dir, '10class', cname, vids(v).name));
                a = load(fullfile(top_dir, 'annot', ...
                    sprintf('%s%.2d.mat', cname, v)));
                mov = box_video(vid, a.annotations, a.width, a.height);
                writer = VideoWriter(fullfile(top_dir, 'combined', ...
                    sprintf('%s%.2d_annot.avi', cname, v)));
                open(writer);
                writeVideo(writer, mov);
            end
        end
    end
end
