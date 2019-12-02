% Source: Matlab user forum discussion
% https://www.mathworks.com/matlabcentral/answers/224604-implementing-gillespie-s-algorithm
% by Harley Day on 11 Feb 2019
h = figure;
hold on;
lineStyles = ['-', '--', ':', '-.'];
markers = ['+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

initPop = [3000 3000; 1000 1000; 100 100; 10 10];
for i=1:length(initPop)
  [t, Y] = eqn38 ( initPop(i, :) );
  c = colors(mod(i-1, length(colors))+1);
  stairs(t, Y, ...
         'LineWidth', 0.4, ...
         'Color', c, ...
         'MarkerSize', 0.5, ...
         'Marker', markers(mod(i-1, length(markers))+1), ...
         'MarkerFaceColor', c, ...
         'MarkerEdgeColor', c);

  legendText{i} = ['$Y_0$ = [' num2str(initPop(i, :)) ']'];
end

line ( [0, 5], [1000, 1000] );
xlim([0 5]);
legendText{length(initPop)+1} = 'Deterministic steady state';
legend( legendText, 'interpreter', 'latex' );

grid on;
hold off;
xlabel('Minutes');
ylabel('Number of Y molecules');
