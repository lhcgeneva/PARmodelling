%% Simulations with real off rates and Ds
sizes = 10:2:14;
s = Simulations({{{'DB', 'DA', 'koffA', 'koffB', 'L', 'StoV'}, {'ratioAB'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {sizes, 15/27})}, {[1.3, 1.4, 1.5]}}}, 1);
s.run_simulations();