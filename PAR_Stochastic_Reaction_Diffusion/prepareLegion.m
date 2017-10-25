c = parcluster('legion_R2016b');
RDSimu = createCommunicatingJob (c, 'Type', 'Pool');
RDSimu.AttachedFiles = {'PARsRDv1_1_StoV.m'
    'align_data.m'
    'create_param_table.m' 
    'assign_slopes.m'       
    'minimize_distance.m'
    'average_datasets.m'
    'plot_Simu.m'
    'Simulations.m'
    'create_change.m'
    'run_simulations.m'};
simulation_duration_ms = 100000;
num_workers = 12;
RDSimu.NumWorkersRange = [num_workers, num_workers];
task = createTask (RDSimu, @runCluster, 0, {});
submit(RDSimu);