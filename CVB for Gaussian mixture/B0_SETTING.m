%///////////////////////////////////////////////////// setting for Ground truth
setting.M = 2;   % number of data's dimension
%-------------------------------------------------
setting.K = 4;   % number of clusters
%-------------------------------------------------
setting.N = 100; % number of time points
%-------------------------------------------------
setting.offset = 1; % offset for origin of x-y plan
%/////////////////////////////////////////////////////

%///////////////////////////////////////////////////// setting for Algorithms
%-------------------------------------------------
setting.maxLoop = 200; % maximum number of iteration
%-------------------------------------------------
setting.init_pos = [-1 0; 0 1; 1 0; 0 -1]'; % matrix of initial position of cluster means (dimension x K)
%-------------------------------------------------
setting.ELBOthresh = 0.01; % converged  if [ELBO < input.ELBOthresh]
%/////////////////////////////////////////////////////

%/////////////////////////////////////////////////////
%-------------------------------------------------
setting.MonteCarlo = 10;    % number of Monte Carlo runs
%-------------------------------------------------
setting.Radius = [0.1,1:8]; % varying Radius (x-axis)
%-------------------------------------------------
setting.plotRadius = 4; % plot the case of Radius = 4
%/////////////////////////////////////////////////////