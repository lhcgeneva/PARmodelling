function plot_Simu( Simu, show_params )
%PLOT_SIMU Plots all Simulations (last timepoint) using tight subplots
%   Also puts the parameters for each performed simulation on the plot that
%   are contained in show_params (e.g. show_params = {'kOn', 'kOff',
%   'kcB'})

numSubPlots = length(Simu.AvSimus);
dimNumR = ceil(sqrt(numSubPlots));
if dimNumR * (dimNumR - 1) >= numSubPlots
    dimNumC = dimNumR -1;
else
    dimNumC = dimNumR;
end

figure;
ts =  tight_subplot(dimNumR, dimNumC, 0.01, 0.01, 0.01);
counter = 1;
for i = 1 : dimNumR
    for j = 1 : dimNumC
        axes(ts(counter));
        hold on;
        plot(Simu.AvSimus{counter}{1}(:, end), 'r');
        plot(Simu.AvSimus{counter}{2}(:, end), 'b');
        %Show legend only in first graph
        if counter == 1
            legend('Species A', 'Species B', 'Location', 'BEST');
        end
        if nargin == 2
            %Get text position, adjust axes
            m = max(max([Simu.AvSimus{counter}{1}(:, end), ...
                    Simu.AvSimus{counter}{2}(:, end)]));%max of A & P
            l = length(Simu.AvSimus{counter}{1}(:, end));%numel in x
            h = m+0.1*m*(length(show_params)+1); %Total height of graph
            axis([0 l 0 h]);
            for i = 1 : length(show_params)
                text(10, h-0.1*m*i, [show_params{i}, ' = ', ...
                    num2str(Simu.param_table_summary.(show_params{i})(counter))]);
                text(10, 10,'Check whether A and B are in correct order!!!');
            end
        end
        counter = counter + 1;
        if counter > numSubPlots; break; end;
    end
end