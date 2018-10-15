classdef Simulations < handle
    %Simulation Contains one set of simulation for specified parameters
    %   Detailed explanation goes here
    
    properties
        Param_list;
        param_table_summary;
        param_table_full;
        num_runs = 10;       %Number of runs per dataset
        keepBinSize=1; %0 If Bin size is fix, 1 if bin size changes depending on system size
        numTpoints;
        SimulationCell;     %Contains simulation structures.
        params;
        simparams;
        AvSimus;
        AllSimus;
        AllSimusLastTimepoint;
        AllSimusAlignedLastT;
        ExperimentArray;
    end
    
    methods
        function Simu = Simulations(Param_list, num_runs)
            if nargin == 2
                Simu.simparams.totalT    = 1000;
                Simu.simparams.boundPos  = 100;
                Simu.simparams.numDom    = 2;
                Simu.simparams.BoundaryC = 'REF';
                Simu.simparams.numTpoints= 100;
                
%                 Assign standard parameters to use as template later on.
                Simu.params.L        = 134.6;
                Simu.params.bins     = 135;
                Simu.params.psi      = 0.174;
                Simu.params.J        = 79;
                Simu.params.Bconc    = 79;
                Simu.params.DA       = 0.28;
                Simu.params.DB       = 0.15;
                Simu.params.alpha    = 1;
                Simu.params.beta     = 2;
                Simu.params.koffA    = 0.0054; % A > cytoA
                Simu.params.konA     = 0.00858; % cytoA > A
                Simu.params.koffB    = 0.0073; % B > cytoB
                Simu.params.konB     = 0.0474; % cytoB > B
                Simu.params.kcA      = 0.19; % A + alpha*B > Acyto + alpha*B
                Simu.params.kcB      = 2; % B + beta*A > Bcyto + beta*A
                Simu.params.p        = 0.0;
                Simu.params.ratioAB  = 1.56;

                %%%%%%%% Symmetric system %%%%%%%%
%                 Simu.params.L        = 134.6/2;
%                 Simu.params.bins     = 135;
%                 Simu.params.psi      = 0.174;
%                 Simu.params.J        = 79;
%                 Simu.params.Bconc    = 79;
%                 Simu.params.DA       = 0.15;
%                 Simu.params.DB       = 0.15;
%                 Simu.params.alpha    = 2;
%                 Simu.params.beta     = 2;
%                 Simu.params.koffA    = 0.0073; 
%                 Simu.params.konA     = 0.0474;
%                 Simu.params.koffB    = 0.0073; 
%                 Simu.params.konB     = 0.0474; 
%                 Simu.params.kcA      = 1; % A + alpha*B > Acyto + alpha*B
%                 Simu.params.kcB      = 1; % B + beta*A > Bcyto + beta*A
%                 Simu.params.p        = 0.0;
%                 Simu.params.ratioAB  = 1;

%                 %%%%%%% Symmetric system paper%%%%%%%%
%                 Simu.params.L        = 120;
%                 Simu.params.bins     = 120;
%                 Simu.params.psi      = 0.3897;
%                 Simu.params.J        = 79;
%                 Simu.params.Bconc    = 79;
%                 Simu.params.DA       = 0.1;
%                 Simu.params.DB       = 0.1;
%                 Simu.params.alpha    = 2;
%                 Simu.params.beta     = 2;
%                 Simu.params.koffA    = 0.005; 
%                 Simu.params.konA     = 0.006;
%                 Simu.params.koffB    = 0.005; 
%                 Simu.params.konB     = 0.006; 
%                 Simu.params.kcA      = 1; % A + alpha*B > Acyto + alpha*B
%                 Simu.params.kcB      = 1; % B + beta*A > Bcyto + beta*A
%                 Simu.params.p        = 0.0;
%                 Simu.params.ratioAB  = 1;
                
                
                Simu.num_runs = num_runs;
                Simu.Param_list = Param_list;
                [Simu.param_table_summary, Simu.param_table_full] = Simu.create_param_table;
%                 Simu.SimulationCell = Simu.run_simulations;
%                 Simu.average_datasets();
%                          
%                 Simu.create_change();
                
                
            end
        end
    end
    
    methods (Access = public)
        run_simulations(Simu);
        [param_table_summary, param_table_full] = create_param_table(Simu);
        average_datasets(Simu);
        assign_slopes(Simu);
        create_change(Simu);
        plot_Simu( Simu, show_params, DetSimus );
    end
    
    methods (Access = public, Static = true)
        [ alignedArray] = align_data(unalignedArray, referenceCurve);
        [ shift, shiftedCurve ] = minimize_distance(curve1, curve2);
    end
end

