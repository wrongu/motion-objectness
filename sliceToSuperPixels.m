function N = sliceToSuperPixels(S, w, h)
    N = zeros(h, w);
    % S(k).points is a 2xP matrix where P is the number of points tracked.
    % Each column is a point [row; col] in the image, i.e. [y; x]
    for sp=1:length(S)
        N(floor(S(sp).points(1,:)+1), floor(S(sp).points(2,:))+1) = sp;
    end
end