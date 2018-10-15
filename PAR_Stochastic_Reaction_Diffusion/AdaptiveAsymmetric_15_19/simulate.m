%% Simulations with real off rates and Ds
sizes = 12:0.5:17;
s = Simulations({{{'DB', 'DA', 'koffA', 'koffB', 'L', 'psi', 'bins'}, {'ratioAB'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {sizes, 15/27}), 2*round(sizes)}, {1.4:0.1:1.9}}}, 1);
s.run_simulations();
%%