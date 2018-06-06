%BURGERSCONTROL demonstrates the control of Burgers PDE using the
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



% the Burgers equation is 
% dv/dt + vdv/dx = nu d^2v/dx^2 
% domain [0,1) with peridoic boundary


% PDE solver parameters
SimPar.N = 100; % spatial grid size
SimPar.dt = 0.01;   % time step
SimPar.N_t = 20;  % length of 
SimPar.nu = 0.01;
x = linspace(0,1,SimPar.N)';    % spatial grid


% some initial conditions 
IC1 = exp(-(((x)-.5)*5).^2); % gaussian initial condition
IC2 = sin(4*pi*x).^2;   % sinusoidal initial condition


% gaussian control profiles
f1 = exp(-(((x)-.25)*15).^2);
f2 = exp(-(((x)-.75)*15).^2);

% bounds on control
umax = .1;
umin =-.1;

%% ************************PART I: Identification ********************** %%
% STEP1: Collecting data for EDMD
Ntraj = 100;    % number of trajectories
SimLength = 200; % Number of steps per trajectory
[X,Y,U] = CollectData(SimPar,IC1,IC2,f1,f2,umin,umax,Ntraj,SimLength);

%% STEP2: EDMD
% choose the observable (to form zeta by delay embedding)
% load BurgersTrajectoryData.mat
CollectMeasurements = @(x) x(10:10:end,:);  
nd = 5; % dimension of embedding (set 1 for no embedding)

% % FULL_STATE observation
% CollectMeasurements = @(x) x(1:end,:);  
% nd = 1; % dimension of embedding (set 1 for no embedding)


[A,B,C,BuildKoopmanState]= SystemID_via_EDMD(X,Y,U,Ntraj,SimLength,CollectMeasurements,nd);

%% ************************PART II: control ********************** %%
% MPC setup for the Koopman linear system
% z'=Az+Bu
% y = Cz
n = size(A,1);  % state dimension
r = size(C,1);  % output dimension
nu= 2;          % input dimension




% cost matrices
Qy = eye(r); % cost matrix on state discrepancy, represents energy
R = eye(nu); % input energy (choosing small values may cause chattering)


% Prediction horizon
Tpred = 0.1;
Np = round(Tpred / SimPar.dt);

% State constraints
xlift_min = nan(n,1);
xlift_max = nan(n,1);

% Input constraints
umin_mpc = umin.*ones(nu,1);
umax_mpc = umax.*ones(nu,1);

% generate the MPC controller handle with no y_ref yet
[~,~,kmpc]  = qpOases_MPC_controller(A,B,C,0,Qy,R,Qy,Np,umin_mpc, umax_mpc, xlift_min, xlift_max,'qpoases',[],[],[],[],[]);

% closed-loop simulation setup
% nonlinear solver
f = @(x,u)(BurgerSolver(x,u(1)*f1+u(2)*f2,SimPar));

Tsim = 6;   % simulation length
Nsim = Tsim / SimPar.dt;
t = 0:SimPar.dt:Tsim;

% initial condition
a = 0.2;
X = a*IC1+(1-a)*IC2;

% reference signal - feel free to play with
xr_t = zeros(size(t));
xr_t(t<=2)=0.5;
xr_t(t>2 & t<4) = 1;
xr_t(t>=4) = 0.5;
Xref = ones(SimPar.N,1)*xr_t;


% initial run to build the initial delay-embeded state
U=[]; 
disp('initialize the closed-loop simulation...')
for i = 1:nd
    X = [ X, f(X(:,end),[0;0]) ];
    U = [U [0;0]]; % choice of input is arbitrary
end

% and now the closed-loop simulation
OutputError=[];
for i = 1:Nsim
    if(mod(i,10) == 0)
            fprintf('Closed-loop simulation, %f %% completed \n', 100 * i / Nsim)
    end

    % reference output
    yref = CollectMeasurements(Xref(:,i));
    
    
    % create the state of the Koopman linear system via embedding and lift
    z = BuildKoopmanState(X,U);

    
    % build the control input
    u = kmpc(z,yref);

    % advance in time and store data
    X = [ X, f(X(:,end),u(:,1)) ];
    U = [U u(:,1)]; 
    
    ynow = CollectMeasurements(X(:,end));
    OutputError = [ OutputError, (ynow-yref)'*Qy*(ynow-yref) ];
end


%% movie presentation

for i = 1:Nsim
    figure(10);
    % state
    subplot(2,1,1)
    plot(x,X(:,i),'color','[0 0 0.85]','linewidth',3); hold on
    plot(x,Xref(:,i),'color','[0.9 0 0]','linestyle','--','linewidth',2.5)
    hold off
    xlabel('$x$','interpreter','latex','fontsize',20);
    ylabel('$v(t,x)$','interpreter','latex','fontsize',20);
    ylim([-0.0,1.2]);

    % input
    subplot(2,1,2)
    plot(x(1:SimPar.N/2),f1(1:SimPar.N/2)*U(1,i),'-r','linewidth',3); hold on
    plot(x(SimPar.N/2:end),f2(SimPar.N/2:end)*U(2,i),'color',[0 0.7 0],'linewidth',3)
    xlabel('$x$','interpreter','latex','fontsize',20);
    ylabel('$u(t,x)$','interpreter','latex','fontsize',20);
    hold off
    ylim([-0.11,0.11]);
    pause(0.001);
end


%% final plots
set(0,'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');
figure(11)
subplot(3,2,1),cla
plot(x,f1); hold on
plot(x,f2); 
xlabel('$x$');
legend({'$f_1$','$f_2$'},'location','North')
title('actuation profiles')


subplot(3,2,2),cla
plot(t, U(1,nd:end)); hold on
plot(t, U(2,nd:end))
xlabel('$t$');
legend({'$u_1$ ','$u_2$'})
title('control inputs')
xlim([0,max(t)])
ylim([min(umin) max(umax)])

subplot(3,2,3),cla
[TT, XX] = meshgrid(t,x);
surf(TT,XX,X(:,nd+1:end)); shading interp
colormap('jet')
xlabel('$t$');
ylabel('$x$');
zlabel('$v$');
title('state evolution')
 view(-18,33)
xlim([0,t(end)])
zlim([0 1.1])


subplot(3,2,4),cla
plot(t,xr_t,'--'), hold on
plot(t,mean(X(:,nd+1:end)));
legend({'reference'})
title('spatial mean')

subplot(3,2,5), hold on
plot(t(2:end),OutputError/size(yref,1),'k--','linewidth',1);
xlabel('$t$');
ylabel('$e$');
xlim([0,t(end)])
title('tracking error')



