% This function draws the bounding box(es) on the given video

function annot_video = box_video(video, annotations, color)
    vprops = get(video, {'Height', 'Width', 'NumberOfFrames'});
    h=vprops{1}; w=vprops{2}; frames=vprops{3};
    assert(frames == length(annotations));
    
    if (nargin < 3)
        color = zeros(1,1,3);
        color(1) = 1;
    else
        color = reshape(color, [1 1 3]);
    end
    
    annot_video(1:frames) = ...
        struct('cdata', zeros(h, w, 3, 'uint8'), 'colormap', []);
    
    for f=1:frames
        % read video frame
        cdata = read(video, f);
        % draw box
        a = annotations{f};
        cdata(a.ytl+1:a.ybr+1, a.xtl, :) = repmat(color, [a.ybr-a.ytl+1 1 1]);
        cdata(a.ytl+1:a.ybr+1, a.xbr, :) = repmat(color, [a.ybr-a.ytl+1 1 1]);
        cdata(a.ytl, a.xtl+1:a.xbr+1, :) = repmat(color, [a.xbr-a.xtl+1 1 1]);
        cdata(a.ybr, a.xtl+1:a.xbr+1, :) = repmat(color, [a.xbr-a.xtl+1 1 1]);
        % save frame
        annot_video(f).cdata = cdata;
    end
end