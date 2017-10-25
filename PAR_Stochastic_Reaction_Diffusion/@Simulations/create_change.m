function  create_change( Simu )
%CREATE_CHANGE Summary of this function goes here
%   Detailed explanation goes here
for i = 1: length(Simu.AvSimus)
    Simu.AllSimusLastTimepoint{i} = cellfun(@(x) x(:, 12, :), ...
                        Simu.AllSimus{i}, 'UniformOutput', false);
    Simu.AllSimusLastTimepoint{i} = cellfun(@(x) squeeze(x), ...
                        Simu.AllSimusLastTimepoint{i}, 'UniformOutput', false);
end
for i = 1 : length(Simu.AvSimus)
    for j = 1:2
        Simu.AllSimusAlignedLastT{i}{j} = Simu.align_data(Simu.AllSimusLastTimepoint{i}{j}, ...
            Simu.AvSimus{i}{j}(:,12));
        [~, dd{i}{j}] = max(smooth(mean(Simu.AllSimusAlignedLastT{i}{j}, 2), 15));
%         side{i}{j}{1} = Simu.AllSimusAlignedLastT{i}{j}(:, dd{i}{j}:dd{i}{j}+29);
            %pPars
        temp = [Simu.AllSimusAlignedLastT{i}{j}; Simu.AllSimusAlignedLastT{i}{j}];
        temp(50+dd{i}{j}-30:50+dd{i}{j})
        side{i}{j}{2} = Simu.AllSimusAlignedLastT{i}{j}(:, dd{i}{j}:dd{i}{j}+29);
    end
end
                
end

