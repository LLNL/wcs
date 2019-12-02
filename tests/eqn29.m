% Source: Matlab user forum discussion
% https://www.mathworks.com/matlabcentral/answers/224604-implementing-gillespie-s-algorithm
% by Harley Day on 11 Feb 2019
function [tplot, Yplot] = eqn29 ( Y0 )
% function to simulate reaction dynamics of equation 29 from Gillespie's
% 1977 paper. The sigle argument to this function specifies the initial
% number of molecules of species Y.
%% Equation 29 reactions
%This gives figure 6 of Gillespie's 1977 paper.
M = 2;                  % Number of reaction pathways
N = 1;                  % Number of molecular species considered
%% STEP 0
% Initial number of molecules of species Y(1) and X(2). In general, a
% vector of length N of initial population numbers for each species.
%------------------------
% equation         | rate
%------------------|-----
% Xbar + Y -> 2*Y  | c1
% Y + Y -> Z       | c2
%------------------------
Y = Y0;
Xbar = 10;
% In general, a vector of length M of stochastic reaction constants.
c = [5/Xbar, 0.005]; % units: (expected number of) reactions per minute
Yplot = Y;            % For plotting.
t = 0;                  % Time variable.
tplot = zeros(1);       % For plotting.
t_stop = 5;
n = 0;                  % Reaction counter.
n_max = 1000;
r_cnt = zeros(M);
%% STEP 1
while t < t_stop
    % Number of molecular reactant combinations available in current state.
    % In general, h is a vector of length M consisting of combinatorial
    % functions of molecular population numbers X (which is of length N).
    h = [Xbar*Y, Y*(Y-1)/2];
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
        case 2
            Y = Y + 1;
        case 1
            Y = Y - 2;
    end
    r_cnt(mu) = r_cnt(mu) + 1;
    n = n + 1;
    % At this point, all the physics has been simulated, it only remains to
    % record the new values of t and X to vectors to we can plot it later.
    Yplot(n+1,:) = Y;
    tplot(n+1) = t;  
end
r_cnt
end
