% This function draws the bounding box(es) on the given video

function annot_video = box_video(video, annotations, a_width, a_height, color)
    vprops = get(video, {'Height', 'Width', 'NumberOfFrames'});
    h=vprops{1}; w=vprops{2}; frames=vprops{3};
    assert(frames == length(annotations));
    
    if (nargin < 5)
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
        cdata = draw_box(cdata, ...
            scale_box(video, [a.xtl a.ytl a.xbr a.ybr], a_width, a_height), color);
        % save frame
        annot_video(f).cdata = cdata;
    end
end

function img = draw_box(img, box, color)
    xtl = box(1);
    ytl = box(2);
    xbr = box(3);
    ybr = box(4);
    img(ytl:ybr, xtl, :) = repmat(color, [ybr-ytl+1, 1, 1]);
    img(ytl:ybr, xbr, :) = repmat(color, [ybr-ytl+1, 1, 1]);
    img(ytl, xtl:xbr, :) = repmat(color, [xbr-xtl+1, 1, 1]);
    img(ybr, xtl:xbr, :) = repmat(color, [xbr-xtl+1, 1, 1]);
end