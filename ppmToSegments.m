% This function takes a 'Segments???.ppm' file output by the Brox/Malik
% motion segmentation function and returns a Kx1 struct, each containing
% a matrix of points for cluster k. The points are 2D column vectors
%
% Author: Richard Lange
% Date Created: March 10, 2013

function Segments = ppmToSegments(ppm_image)
% if(nargin < 2)
%     f_in = figure();
% end
I = ppm_image;
[h w ~] = size(I); % height, width, and image depth
I_accounted = zeros(h,w);
Segments = struct('points', {[], [], []});
next_k = 1;
% a dictionary-like sparse matrix for looking up k from color
k_lookup = sparse(256^3,1);
r = 1;
c = 1;
for r=1:h
    for c=1:w
        if(~I_accounted(r,c))
            val = double(I(r,c,:));
            if any(val ~= 255)
                lookup_index = (val(1)-1)*256^2 + (val(2)-1)*256+val(3) + 1;
                k = k_lookup(lookup_index,1) + 0; % the plus 0 converts from a sparse entry to a scalar value
                if k == 0
                    k = next_k;
                    fprintf('new cluster: %d\n', k);
                    k_lookup(lookup_index,1) = next_k;
                    next_k = next_k+1;
                end
                % find the center associated with this non-white pixel:
                [r_c c_c] = get_center(I, I_accounted, r, c);
                Segments(k).points = horzcat(Segments(k).points, [r_c; c_c]);
                I_accounted( max(r_c-4,1) : min(r_c+4, h), ...
                    max(c_c-4,1) : min(c_c+4, w) ) = 1;
%                 fprintf('found center at %d, %d, cluster %d\n', r_c, c_c, k);
%                 figure(f_in);
%                 imshow(I_accounted);
%                 pause;
            end
        end
    end
end
disp(k_lookup);
end

function [r_ c_] = get_center(I, I_acc, r, c)
for r_ = max(r-2,1):min(r+2,size(I,1))
    for c_ = max(c-2,1):min(c+2,size(I,2))
        if ~I_acc(r_, c_) && is_center(I, r_, c_)
            return;
        end
    end
end
end

function c = is_center(I, r, c)
% check that full block centered at r,c is set:
c = all( ...
    all( ...
    any( ...
    I( max(r-2,1):min(r+2,size(I,1)),  max(c-2,1):min(c+2,size(I,2)), : ) ~= 255, 3)));
end