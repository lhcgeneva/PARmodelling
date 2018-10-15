%% Simulations with real off rates and Ds
sizes = 12:0.5:15;
s = Simulations({{{'DB', 'DA', 'koffA', 'koffB', 'L', 'psi', 'bins'}, {'ratioAB'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {sizes, 15/27}), 2*round(sizes)}, {1.5:0.1:1.9}}}, 10);
s.run_simulations();
save s1.mat
%%
%%
a = s.AvSimus{1}{1}(:, end);
mean(a(1:length(a)/2))