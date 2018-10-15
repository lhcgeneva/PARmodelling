function [ param_table_summary, param_table_full ] = create_param_table(Simu)
%create_param_table Takes list of parameters and creates parameter table
%   Each row in the output table corresponds to a set of parameters for
%   which the simulation is run. The parameter combinations are constructed
%   by first creating two arrays (N_R and N_B) which contain, in the order
%   of input, the number of times each parameter has to be repeated in
%   order to be able to accomodate all the still following parameters, as well
%   as how many times each of these has to be written underneath each
%   other, in order to fill all the ones that have already been filled.
%   Example: input: {{'L', 'psi'}, {[1, 2, 3], [0.1, 0.2]}}.
%   table is supposed to look like this:
%   L       psi     ...
%   1       0.1     ...                     Block 1 for psi
%   1       0.2     ...  Block 1 for L      Block 2 for psi
%   ----------------------------            
%   2       0.1     ...                     Block 1 for psi
%   2       0.2     ...  Block 2 for L          .
%   ----------------------------                .
%   3       0.1     ...                         .
%   3       0.2     ...  Block 3 for L
%   ----------------------------
%   N_B (number of repetitions per block, for L this would be 2, for psi 1)
%   N_R (number of times all blocks of one parameter are repeated, this
%   would be 1 for L and 3 for psi) 
%   Both of these arrays are created by calculating all the
%   combinations possible for the parameters that have yet to be
%   incorporated (N_B) and the ones that have already been incorporated
%   (N_R) N_B is calculated cumprod(n1, n2, n3 ...) were n is the length of each
%   parameter's input array (3 for L, 2 for psi). N_R is calculated as the
%   'inverse' cumulative product. N_B has to be shifted by 1 (see for loop)
%   because the number of upcoming combinations is important, not including
%   the current parameter's array length.

% param_table(prod(cellfun(@length, Simu.Param_list{2})), :) = struct2table(Simu.params);
param_table(prod(cellfun(@length, cellfun(@(x) x{1}, Simu.Param_list{2}, ...
            'UniformOutput', false))), :) = struct2table(Simu.params);
N_R = [1, cumprod(cellfun(@length, cellfun(@(x) x{1}, Simu.Param_list{2}, ...
            'UniformOutput', false)))];
N_B = fliplr(cumprod(fliplr([cellfun(@length, cellfun(@(x) x{1}, Simu.Param_list{2}, ...
            'UniformOutput', false)), 1])));

for j = 1 : N_B(1)
    param_table(j, :) = struct2table(Simu.params);
end

for p = 1 : length(Simu.Param_list{1})
temp_index = 1;
    for i = 1 : N_R(p)
        for k = 1:length(Simu.Param_list{2}{p}{1})
            for j = 1:N_B(p+1)
                for l = 1 : length(Simu.Param_list{1}{p})
                    param_table.(Simu.Param_list{1}{p}{l})(temp_index) = Simu.Param_list{2}{p}{l}(k);
                end
                temp_index = temp_index + 1;
            end
        end
    end
end

% Add column with linear index to front.
LinInd = (1:N_B(1))';
LinInd = table(LinInd);
param_table_summary = [LinInd, param_table];
param_table_full = repmat(param_table_summary, Simu.num_runs, 1);
param_table_full.LinInd = (1:(Simu.num_runs*N_B(1)))';
end

