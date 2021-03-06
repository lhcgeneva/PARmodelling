s = Simulations({{{'DB', 'DA'}}, {{linspace(6,20,12), linspace(6,20,12)}}}, 5);
s.run_simulations;
save diffusion_6_20;
s = Simulations({{{'DB', 'DA'}}, {{linspace(0.005,0.09,12), linspace(6,20,12)}}}, 5);
s.run_simulations;
save diffusion_0005_009;

%% Simulations with real off rates and Ds
sizes = 100:10:200;
s = Simulations({{{'DB', 'DA', 'koffA', 'koffB', 'L', 'StoV'}},...
                {{fDP2(sizes), fDP6(sizes), fOP2(sizes), fOP6(sizes),...
                  sizes,S_to_V('Circumference', {sizes, 15/27})}}}, 1);
s.run_simulations();