s = Simulations({{{'DB', 'DA'}}, {{linspace(0.1,5,12), linspace(0.1,5,12)}}}, 5);
s = Simulations({{{'DB', 'DA'}}, {{linspace(6,20,12), linspace(6,20,12)}}}, 5);
s = Simulations({{{'koffB', 'koffA'}}, {{linspace(6,20,12), linspace(6,20,12)}}}, 5);
s.run_simulations();
save('/home/ucgahub/Scratch/Matlab_remote_jobs/testRD.mat');