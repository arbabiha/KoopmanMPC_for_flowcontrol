function [X,Y,U] = CollectData(SimPar,IC1,IC2,f1,f2,umin,umax,Ntraj,SimLength)
% Collect data for Burgers equation
disp('Starting data collection ...')


% random input
Ubig= rand(2,SimLength,Ntraj)* (umax-umin)+umin;

% Transition mapping of the controleld dynamical system
% control is a linear combination of u1 and u2
f = @(x,u)(BurgerSolver(x,u(1)*f1+u(2)*f2,SimPar));


% initialize 
X = []; Y = []; U=[];

% loop pver trajectories
for i = 1:Ntraj
    xx = [];
    % Intial state is a random convex combination of IniitialConditions
    b = rand;
    a = [b,1-b];
    xx =b*IC1 + (1-b)*IC2;
    tic
    fprintf('Trajectory %d out of %d \n',i,Ntraj)
    % loop over each time step
    for j = 1:SimLength 
        xx = [xx f(xx(:,end),Ubig(:,j,i))];
        U  = [U,Ubig(:,j,i)];
        % if the solution diverges, go to the next trajectory
        if ~isempty(find(isnan(xx(:,end)),1))
            break
        end
        
    end
    toc
    % Store
    X = [X xx(:,1:end-1)];
    Y = [Y xx(:,2:end)];

end


save('BurgersTrajectoryData','X','Y','U','SimLength','Ntraj')
end