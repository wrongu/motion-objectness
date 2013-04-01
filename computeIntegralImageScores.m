function score = computeIntegralImageScores(integralImage,windows)

windows = round(windows);
%windows = [xmin ymin xmax ymax]
%computes the score of the windows wrt the integralImage
height = size(integralImage,1);
index1 = height*windows(:,3) + (windows(:,4) + 1); % bottom-right
index2 = height*(windows(:,1) - 1) + windows(:,2); % top-left
index3 = height*(windows(:,1) - 1) + (windows(:,4) + 1); % bottom-left
index4 = height*windows(:,3) + windows(:,2); % top-right
% score = BR + TL - BL - TR
%   gives the sum of the values in the integral image _within_ the windows
score = integralImage(index1) + integralImage(index2) - integralImage(index3) - integralImage(index4);

end