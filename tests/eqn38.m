% Source: Matlab user forum discussion
% https://www.mathworks.com/matlabcentral/answers/224604-implementing-gillespie-s-algorithm
% by Harley Day on 11 Feb 2019

%% Lotka reactions
%--------------------------
% equation           | rate
%--------------------|-----
% Xbar + Y1 -> 2*Y1  |  c1
% Y1 + Y2 -> 2*Y2    |  c2
% Y2 -> Z            |  c3
%--------------------------
function [tplot, Y1plot, Y2plot] = eqn38 ( Y0 )
M = 3;                  % Number of reaction pathways
N = 2;                  % Number of molecular species considered

%% STEP 0
% Initial number of molecules of species Y1 and Y2.
%------------------------
Y1 = Y0(1);
Y2 = Y0(2);
Xbar = 10;

% In general, a vector of length M of stochastic reaction constants.
c = [10/Xbar, 0.01, 10]; % units: (expected number of) reactions per minute
Y1plot = Y1;             % For plotting.
Y2plot = Y2;
t = 0;                   % Time variable.
t_max = 30;
tplot = zeros(1);        % For plotting.
n = 0;                   % Reaction counter.
n_max = 1000000;
r_cnt = zeros(M);

%% STEP 1
while t < t_max
  % Number of molecular reactant combinations available in current state.
  % In general, h is a vector of length M consisting of combinatorial
  % functions of molecular population numbers Xbar (which is of length N).
  h = [Xbar*Y1, Y1*Y2, Y2];
  % a is the propensity of the reaction pathway in current state. In
  % general, a vector of length M
  a = h.*c;
  % a0 is total propensity that anything happens. This number emerges
  % more out of mathematical necessity than physical intuition.
  a0 = sum(a);

  %% STEP 2
  r = rand(2,1);
  tau = -log(r(1))/a0;
  % Note 0<=mu<=M. This decides which reaction occurs. If mu=0, no
  % reaction occurs. The switch statement does nothing in this case.
  mu = sum(r(2)*a0 <= cumsum(a));
    
  %% STEP 3
  t = t + tau;
  % Adjust population levels based on reaction formula(s).
  switch mu
    case 3
        Y1 = Y1 + 1;
    case 2
        if (Y1 > 0)
            Y1 = Y1 - 1;
            Y2 = Y2 + 1;
        else
            disp('reaction 2 cannot fire');
        end
    case 1
        if (Y2 > 0)
            Y2 = Y2 - 1;
        else
            disp('reaction 3 cannot fire');
        end
  end
  r_cnt(mu) = r_cnt(mu)+1;
  n = n + 1;
  % At this point, all the physics has been simulated, it only remains to
  % record the new values of t and Xbar to vectors to we can plot it later.
  Y1plot(n+1) = Y1;
  Y2plot(n+1) = Y2;
  tplot(n+1) = t;
end
r_cnt
end
