function [ shift, shiftedCurve ] = minimize_distance( curve1, curve2 )
%MINIMIZE_DISTANCE ...between curve1 and curve2. Returns shift + shifted curve 
%   Detailed explanation goes here

dist = @(x_off) sum(([curve1(end-x_off+1:end);curve1(1:end-x_off)]-curve2).^2);
%For-loop because arrayfun is slower
for i = 1 : length(curve1)
    distArray(i) = dist(i);
end
[~, shift] = min(distArray);
shiftedCurve = [curve1(end-shift+1:end);curve1(1:end-shift)];
end

