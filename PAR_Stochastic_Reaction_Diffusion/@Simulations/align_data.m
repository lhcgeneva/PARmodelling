function [ alignedArray ] = align_data( unalignedArray, referenceCurve )
%ALIGN_DATA Summary of this function goes here
%   Detailed explanation goes here
alignedArray = zeros(size(unalignedArray));
[~, num_sets] = size(alignedArray);
for i = 1:num_sets
    [~, alignedArray(:, i)] = Simulations.minimize_distance(unalignedArray(:, i) ,...
                                            referenceCurve);
end

end

