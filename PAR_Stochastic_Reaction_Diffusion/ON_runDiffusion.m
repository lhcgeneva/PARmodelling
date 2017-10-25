s = Simulations({{{'DB', 'DA'}}, {{linspace(6,20,12), linspace(6,20,12)}}}, 5);
s.run_simulations;
save diffusion_6_20;
s = Simulations({{{'DB', 'DA'}}, {{linspace(0.005,0.09,12), linspace(6,20,12)}}}, 5);
s.run_simulations;
save diffusion_0005_009;