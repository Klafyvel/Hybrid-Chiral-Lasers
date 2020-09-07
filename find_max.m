% find_max.m
%
% Find the local maxima of a 2D array and return their row and column. A local maxima is identified as being greater than all values around it. 
% Values outside the scope are considered to be infinite.

function [row,col] = find_max(array, radx, rady)
    [w,h] = size(array);
    right = circshift(array, [-1,0]);
    right(:,end) = NaN(w,1);
    left = circshift(array, [1,0]);
    left(:,1) = NaN(w,1);
    above = circshift(array, [0,1]);
    above(1,:) = NaN(1, h);
    below = circshift(array, [0,-1]);
    below(end,:) = NaN(1, h);
    [row_t, col_t] = find((array > left) & (array > right) & (array > above) & (array > below));

    % In order to detect peaks, only take maximas that are maximum for the area surrounding.
    k = 1;
    [max_row, max_col] = size(array);
    row = [];
    col = [];
    for j = 1:length(row_t)
        r = row_t(j);
        c = col_t(j);
        max_area = array(max(1,r-radx):min(max_row,r+radx),max(1,c-rady):min(max_col,c+rady));
        if array(r,c) >= max_area
            col(k) = c;
            row(k) = r;
            k = k+1;
        end
    end
end

