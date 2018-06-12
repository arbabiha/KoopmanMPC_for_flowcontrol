%CAVITYCONTROL demonstrates the control of lid-driven cavity flow using the
%Koopman-MPC framework described in 
% "A data-driven Koopman model predictive control framework for nonlinear
% flows" H. Arbabi, M. Korda and I. Mezic

clc,clear
addpath('./thehood')
% quadratic programming solver
addpath('./thehood/qpOASES-3.1.0/interfaces/matlab')


if (exist('qpOASES_sequence','file') ~= 3)    

    error(['You have to activate the MATLAB interface for qpOASES first:' ...
        ' unzip the qpOASES-3.1.0 in "thehood" folder, then '...
        'go to ".\thehood\qpOASES-3.1.0\interfaces\matlab" and run make.m,' ...
        'then run again'])

end



% simulation and flow parameters
SimPar.Reynolds = 13000;
SimPar.dt =.01;            % time step of the ODE solver
SimPar.N = 49;             % size of the computaional grid is (N+1)^2
SimPar.T = 0.2;            % length of each run - continuous-time cavity flow is approximated 
                           % using T-time discrete map




load('CavityStateLibrary.mat','LimitCycle_Re13k','FixedPoint_Re10k','UnstableFixedPoint_Re13k')

% initial condition
 x0 = LimitCycle_Re13k(:,25);   % some point on the limit cyle




%% ************************PART I: Identification ********************** %%
%% STEP2: Identification via EDMD
% there are two options:
% 1- generate the data via runing the following program
% (for default values takes ~10hrs on a powerful desktop)
% DataFileName=GenerateCavityData;

% 2- download the data file cited in README.md and set
DataFileName = 'Cavity_data_4EDMD_0.mat';

% Generating Koopman-linear predictors with various number of measurements
CavitySystemID(DataFileName);


%  load any of the predictors
load('KoopmanLinSys_Re13_k50_random');    % predictor with k=50 random measurements

%% ************************PART II: control ********************** %%
% set up MPC controller
disp('setting up MPC controller ...')

n = size(A,1);  % state dimension
r = size(C,1);  % output dimension


% cost matrices
Qy = speye(r);  % output weight matrix
R = 0; % zero weight on input (there are input constraints)


% Prediction horizon
Tpred = 10;
Np = round(Tpred / SimPar.dt);

% State constraints
x_min = nan(n,1);
x_max = nan(n,1);


% Input constraints
umin = 11/13;
umax = 15/13;


% Precompute and possibly save MPC big matrices to speed construction of the controller. 
% (This "lifting" has nothing to do with the koopman lifting)
[Ab, Bb] = createMPCmatrices(A,{B},Np); 

% compute controller handle
tic
[~,~,mpcCont_lift]  = qpOases_MPC_controller(A,B,C,0,Qy,R,Qy,Np,umin, umax, x_min, x_max,'qpoases',[],[],[],[],[],Ab,Bb);
toc


% reference state
xref = FixedPoint_Re10k;    % unstable fixed point at Re=13k
xref_mpc = xref - x_mean; % mean subtracted! 

% closed-loop simulation setup
% nonlinear solver
f = @(x,u)(NonlinearFlowSolver(x,CreateLidVelocity(u,SimPar.N),SimPar));



%% initialize
X = x0;
U = [];
cost = (CollectOutput(x0 - xref))'*Qy*(CollectOutput(x0 - xref));
Q_KE = getCostMatrix(SimPar.N); % the weight matrix to compute kinetic energy
KE_discrepancy = (x0 - xref)'*Q_KE*(x0 - xref);



% run to build the initial delay embeded state
for i = 1:nd

%     Nonlinear simulation with base input
    X = [ X, f(X(:,end),1) ];
    U = [U 1]; 
end

% now build the initial embedded state
disp('initialization complete')


%% closed-loop simulation

Tsim = 100;
Nsim = Tsim / SimPar.T;

for i = 1:Nsim
    if mod(i,50)==0
        fprintf('Closed-loop simulation, %f %% completed \n', 100 * i / Nsim)
    end
    % create reference output
    yref = CollectOutput(xref_mpc(:,end));

    
    % build the state of the Koopman predictor
    y = CollectOutput(  bsxfun(@minus,X(:,end-nd+1:end),x_mean) ); % collect the measurements 
    y = DelayEmbed(y,nd);   % delay embed the measurements
    ue = DelayEmbed(U(:,end-nd+1:end-1),nd-1);  % delay embed the input
    znow = [y;ue;KE_embed(y);1];  % form the state (g(zeta) in the paper)

    
    % compute control
    u = mpcCont_lift(znow,yref);

    % Nonlinear simulation
    X = [ X, f(X(:,end),u(:,1)) ];
    U = [U u(:,1)];
    
    % compute the tracking error and kinetic energy of state discrepancy
    ynow = CollectOutput(  bsxfun(@minus,X(:,end),x_mean) );
    cost = [ cost, (ynow-yref)'*Qy*(ynow-yref) ];
    KE_discrepancy = [KE_discrepancy, (X(:,end) - xref)'*Q_KE*(X(:,end) - xref)];

end


%% plots
set(0,'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');
t = (0:Nsim)*SimPar.T;
figure(20),clf
subplot(2,2,1)
plot(t,U(1,nd:end),'linewidth',2);
xlabel('$t$','fontsize',12);
title('input')
ylim([umin umax])
set(gca,'YTick',linspace(umin,umax,5))

subplot(2,2,2)
plot(t,KE_discrepancy); hold on
title('kinetic energy of discrepency')

subplot(2,2,3)
plot(t,cost); hold on
title('control tracking error')

subplot(2,2,4)
PlotVorticity(X(:,end)-xref);
hold on
[ Grid ] = CavityGridOperators( SimPar.N );
plot(CollectOutput(Grid.xx),CollectOutput(Grid.yy),'x')
axis square
axis([-1 1 -1 1])
title('vorticity discrepancy at final time and sensor locations')