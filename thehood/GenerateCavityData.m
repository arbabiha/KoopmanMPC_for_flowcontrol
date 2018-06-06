function DataFileName=GenerateCavityData();
% generating trajectories of cavity flow with random initial condition 
% and sequence of random inputs

clear NonlinearFlowSolver CreateLidVelocity




% simulation and flow parameters 
SimPar.Reynolds = 13000; % base Reynolds
SimPar.N = 49;            % size of the computaional grid is (N+1)^2

% bounds on control input
umax = 10/13;
umin = 16/13;

% Time stepping
SimPar.dt=.01;            % time step of the ODE solver
SimPar.np= 100;            % every np*dt seconds the plot is updated

% T is the time increment, it should be a multiple of dt
SimPar.T = 0.2;           
% data collection parameters
Ntraj = 500;
SimLength = 200; % Number of steps per trajectory





load('CavityStateLibrary.mat','LimitCycle_Re13k','FixedPoint_Re10k')

%% ******************************* EDMD ***********************************
% number of inputs
nu = 1;

% Transition mapping of the controleld dynamical system
f = @(x,u)(NonlinearFlowSolver(x,CreateLidVelocity(u,SimPar.N),SimPar)); % 1 input



tic
disp('Starting data collection for cavity flow')

SimLengthVar = SimLength*ones(1,Ntraj);



% Generate control inputs
ContPar.umax = umax;
ContPar.umin = umin;
Ubig = rand(1,SimLength*Ntraj)*(umax-umin) + umin;
Ubig = reshape(Ubig(:,1:SimLength*Ntraj),[1,SimLength,Ntraj]);



%% simulation and data collection


X = []; Y = []; U = [];

for i = 1:Ntraj
    xx = [];
    b = rand;
    % Intial state
    a = [b,1-b]+0.1*randn(1,2);
    xx =a(1)*LimitCycle_Re13k(:,randi(50)) + a(2)*FixedPoint_Re10k;



    if(SimPar.Reynolds <= 300)
        xx = zeros((SimPar.N+1)^2,1); % TO BE REPLCED WITH A LINEAR COMBINATION OF TWO MEAN FLOWS
    end
    
    
    % Simulate one trajectory
    tic
    fprintf('Trajectory %d out of %d \n',i,Ntraj)
    for j = 1:SimLength
        %         fprintf('Trajectory %d out of %d, step %d out of %d \n',i,Ntraj,j,SimLength)
        xx = [xx f(xx(:,end),Ubig(:,j,i))];
        U  = [U,Ubig(:,j,i)];
        % if the solution diverges ...
        if ~isempty(find(isnan(xx(:,end)),1))
            SimLengthVar(i)=j;
            break
        end
        
    end
    toc
    % Store
    X = [X xx(:,1:end-1)];
    Y = [Y xx(:,2:end)];
    
end

DataFileName='Cavity_data_4EDMD.mat';
save(DataFileName,'X','Y','U','SimPar','Ntraj','SimLengthVar','-v7.3')

end
