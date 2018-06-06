function [FinalState,T_pass,KE] = NonlinearFlowSolver(STATE,INPUT,SimulationParameters)
% NonlinearFlow solves the lid-driven cavity flow
% Inputs:
%     STATE: Vector of stream function values on 2D chebyshev points
%     INPUT: Vector of velocity profile on 1D chebyshev points on the top
%     lid
%     SimulationParameters: A struct with the following fields:
%                           Reynolds: Reynolds number , 
%                           N: Order of Chebyshev polynomials in each direction
%                           dt: Time step for simulation (optional)
%                           T: simulation time span



% Outputs:
       
%   FINALSTATE: Vector of stream function values after time T
%   T: elapsed time







    
% accessing auxiliary files and clear persistent memory
clear RHS_cavity_controlled


% blurt out
% disp(['cavity flow at Re=',num2str(SimulationParameters.Reynolds)])

  


% variables that we save in the function memory for repeated use
persistent Grid Operators_q Operators_psi %LHS

if isempty(Grid)
    [Grid]=CollocationGrid_q(SimulationParameters.N);    % the computational grid
    [~,Grid.wc] = clencurt(SimulationParameters.N);Grid.W = kron(Grid.wc,Grid.wc);   % integration weights
     disp('grid created')
end

if isempty(Operators_q)
    [ Operators_q ] = CreateOperators_q( Grid,eye(SimulationParameters.N+1) );
     disp('operators for q created')
end
if isempty(Operators_psi)
    [ Operators_psi ] = CreateOperators_psi( Grid.D,eye(SimulationParameters.N+1));
    disp('operators for psi created')
end

% We normalzie the velocity field using max-norm and pass its max to the
% Reynolds number
NominalInput= (1-Grid.x.^2).^2;
VelocityRatio = max(INPUT)/max(NominalInput);
INPUT = INPUT./VelocityRatio;

% scale the Reynolds number
SimulationParameters.Re=SimulationParameters.Reynolds*VelocityRatio/2;

% if isempty(LHS)
    LHS = Operators_q.Del2 - (0.5*SimulationParameters.dt / SimulationParameters.Re) * Operators_q.Del4;
    LHS(Grid.Walls,:)=0; LHS(Grid.Walls,Grid.Walls)=eye(4*SimulationParameters.N);  % Dirichlet BC on q(x,y)
    LHS = inv(LHS);
% end




% time stepping parameters
  nt = floor(SimulationParameters.T/SimulationParameters.dt);
  Time = (1:nt)*SimulationParameters.dt;
  KE   = zeros(nt,1);
%   PauseStep = SimulationParameters.np;




% initialize
Psi.Previous = STATE;
Psi.Current =  STATE;
Psi.Temp   = zeros(size(STATE));



% disp('iterating ...')

for it=1:nt
    
            % solving the system
            RHS = RHS_cavity_controlled( Operators_psi,Grid,Psi,SimulationParameters,INPUT);
            q_psi = LHS*RHS;
            Psi.Temp = ((1-Grid.xx.^2).*(1-Grid.yy.^2)).*q_psi;
            
            % shifting
            Psi.Previous = Psi.Current;
            Psi.Current  = Psi.Temp;
            


end

% assigning the output
% KE(end) = CalculateKineticEnergy(Operators_psi,Grid,Psi.Current);
FinalState = Psi.Current;
T_pass = Time(it);


end


function [ Operators ] = CreateOperators_q( Grid,I)
%CREATEOPERATORS_q creates the Operator acting on q(x,y)


Dc2 = (kron(Grid.D^2,I)*diag(1-Grid.xx.^2))*(kron(I,Grid.D^2)*diag(1-Grid.yy.^2));
Dx4 = kron(Grid.D4,I)*diag(1-Grid.yy.^2); Dy4 = kron(I,Grid.D4)*diag(1-Grid.xx.^2);
Dxx = kron(Grid.D^2,I)*diag(1-Grid.xx.^2);
Dyy = kron(I,Grid.D^2)*diag(1-Grid.yy.^2);

Operators.Del2 = diag(1-Grid.xx.^2)*Dyy + diag(1-Grid.yy.^2)* Dxx;

Operators.Del4 = Dx4  +   Dy4     +   2*Dc2;

S = (1-Grid.yy.^2).*(1-Grid.xx.^2);


Operators.Dx= spdiags( (-2)*Grid.xx.*(1-Grid.yy.^2),0,length(Grid.xx),length(Grid.xx)) ...
             +bsxfun(@times,S,kron(Grid.D,I));
         
Operators.Dy= spdiags( (-2)*Grid.yy.*(1-Grid.xx.^2),0,length(Grid.xx),length(Grid.xx)) ...
             +bsxfun(@times,S,kron(I,Grid.D));        

end


function [KineticEnergy] = CalculateKineticEnergy(Operators_psi,Grid,StreamFunction)


                ux =  Operators_psi.Dy * StreamFunction;  % u-velocity
                uy =  Operators_psi.Dx * StreamFunction;  % v-velocity: minus removed to reduce simulation time
                KineticEnergy = 0.5 * (Grid.W* (ux.^2) ) + 0.5 * (Grid.W*(uy.^2));
                
end