% Convert video (avi) to directory with frames as ppm images and bmf file
% (see specification for motion segmentation in Brox algorithm)
% returns path of the directory

function dest = videoToBMF(video_path, name, dest)
    dest = fullfile(dest, name);
    if ~exist(dest, 'dir')
        mkdir(dest);
    end
    
    V = VideoReader(video_path);
    vData = get(V, {'NumberOfFrames'});
    num_frames = vData{1};
    
    Vmatrix = read(V);
    
    bmf_text = sprintf('%d %d\n', num_frames, 1);
    
    for f = 1:num_frames
        img = Vmatrix(:,:,:,f);
        img_name = sprintf('%s_%.3d.ppm', name, f);
        fprintf('writing img to %s\n', fullfile(dest, img_name));
        imwrite(img, fullfile(dest, img_name), 'PPM');
        bmf_text = sprintf('%s%s\n', bmf_text, img_name);
    end
    
    bmf_file = fopen(fullfile(dest, sprintf('%s.bmf', name)), 'w+');
    fprintf(bmf_file, bmf_text);
    fclose(bmf_file);
end