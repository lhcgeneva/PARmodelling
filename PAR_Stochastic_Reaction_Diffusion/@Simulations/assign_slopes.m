function assign_slopes( Simu )
%ASSIGN_SLOPE Assigns slope to every Structure A and B in SimulationCell
%   Uses calc_slope 
    for i = 1:numel(Simu.SimulationCell)
        windowSize = 18.1 / (Simu.SimulationCell{i}.params.L / ...
                            Simu.SimulationCell{i}.params.bins);
        windowSize = 2 * round(windowSize/2) + 1;
        Simu.SimulationCell{i}.slopeA = ...
                        calc_slope(flipud(Simu.SimulationCell{i}.A(:, end)), windowSize, 5, 0);
        Simu.SimulationCell{i}.slopeB = ...
                        calc_slope(Simu.SimulationCell{i}.B(:, end), windowSize, 5, 0);
    end
end

