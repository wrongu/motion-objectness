% return box coordinates that fit on the video V, taken from annotations
% that have width a_width and height a_height (w and h from 'annotation')
% box and a_box should be [xtl ytl xbr ybr]
function box = scale_box(V, a_box, a_width, a_height)
    vdata = get(V, {'Width', 'Height'});
    w = vdata{1};
    h = vdata{2};
    box = a_box;
    box(:, [1 3]) = box([1 3]) * w / a_width;
    box(:, [2 4]) = box([2 4]) * h / a_height;
end