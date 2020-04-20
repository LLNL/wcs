% Source: Matlab user forum discussion
% https://www.mathworks.com/matlabcentral/answers/224604-implementing-gillespie-s-algorithm
% by Harley Day on 11 Feb 2019
lineStyles = ['-', '--', ':', '-.'];
markers = ['+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

sample_interval = 0.1; % Sampling interval used
tmax = 5; % Simulation end time used
Y0 = [10 1000 3000]; % Inital Y population used
ssa = 0; % SSA method used: 0 for direct 1 for NRM
seeds = [59 67 73 83 97 103 109 127 137 149 157];
path = ['result-sample-' num2str(sample_interval) '/toplot/'];

FigW=10;
FigH=8;

for s=1:length(seeds)
  seed = seeds(s);
  close all;
  h1 = figure;
  hold on;
  for i_Y0 = 1:length(Y0)
    filename = [ path 'eq29.Y-' num2str(Y0(i_Y0)) '.m-' num2str(ssa) '.seed-' num2str(seed) ];
    cms_idx = i_Y0; % colori-marker-style index
    
    tY = load([ filename '.txt' ]);
    t = tY(:,1);
    Y = tY(:,2);
    clear tY;
    
%    clf;
%    set(h1,'defaulttextinterpreter','latex',...
%                'PaperUnits','centimeters','PaperSize',[FigW FigH],...
%                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
%                'Position',[1,01,FigW,FigH]);

    c1 = colors(mod(cms_idx-1, length(colors))+1);
    h1 = plot(t, Y, ...
             'LineWidth', 0.4, ...
             'Color', c1, ...
             'MarkerSize', 0.5, ...
             'Marker', markers(mod(cms_idx-1, length(markers))+1), ...
             'MarkerFaceColor', c1, ...
             'MarkerEdgeColor', c1);
    clear t;
    clear Y;
    lgn1{i_Y0} = ['$Y_0$ = [' num2str(Y0(i_Y0)) ']'];
  end

  grid on;
  hold off;
  xlim([0 tmax]);
  xlabel('Minutes');
  %ylabel('Number of Y molecules', 'interpreter', 'latex');
  ylabel('Number of Y molecules');

  legend(lgn1, 'interpreter', 'latex', 'Location', 'NorthEast');

  print([filename '-t-Y.pdf'], '-dpdf');
end
