function average_datasets( Simu )
%AVERAGE_DATASETS takes all the created datasets and averages...
%   ...the ones using the same input parameters.
h = height(Simu.param_table_summary);
Simu.AvSimus = cell(h, 1);
Simu.AllSimus = cell(h, 1);

for i = 1 : h 
    dummy = [Simu.SimulationCell{i, :}];
    sz = size([dummy(1).A]);
    Simu.AllSimus{i}{1} = reshape([dummy.A], sz(1), ...
                    sz(2), Simu.num_runs);
    Simu.AvSimus{i}{1} = mean(Simu.AllSimus{i}{1}, 3);
    dummy = [Simu.SimulationCell{i, :}];
    Simu.AllSimus{i}{2} = reshape([dummy.B], sz(1), ...
                    sz(2), Simu.num_runs);
    Simu.AvSimus{i}{2} = mean(Simu.AllSimus{i}{2}, 3);
end

