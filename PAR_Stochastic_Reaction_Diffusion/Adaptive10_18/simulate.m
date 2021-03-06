%% Simulations with real off rates and Ds
sizes = 10:2:18;
s = Simulations({{{'DB', 'DA', 'koffA', 'koffB', 'L', 'StoV'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {sizes, 15/27})}}}, 1);
s.run_simulations();