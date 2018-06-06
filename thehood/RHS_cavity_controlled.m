  function [ RHS ] = RHS_cavity_controlled( Operators_psi,Grid,Psi,SimPar,U_lid )
%RHS_CAVITY Summary computes the RHS for cavity problem

% Advection terms (explicit in time), Diffusion (implicit) and time
% increment
persistent w0 w1 u0 u1 v0 v1

if isempty(w0)
    % On the first iteration, w0 will be empty
%     disp('creating persistent variables in RHS')
    w0 = Operators_psi.del2 * Psi.Previous;      % minus vorticity
    u0 = Operators_psi.Dy *Psi.Previous;
    v0 = (-1)*(Operators_psi.Dx * Psi.Previous);
else
    w0 = w1;
    u0 = u1;
    v0 = v1;
end

    w1 = Operators_psi.del2 * Psi.Current;
    u1 = Operators_psi.Dy *Psi.Current;
    v1 = (-1)*(Operators_psi.Dx * Psi.Current);  % the minus sign is taking time?




RHS = ((-1.5)*SimPar.dt)*(u1.*(Operators_psi.Dx * w1) + v1.*(Operators_psi.Dy * w1) ) +...
      +(0.5*SimPar.dt)*(u0.*(Operators_psi.Dx * w0) + v0.*(Operators_psi.Dy * w0) ) +...
      +(0.5*SimPar.dt/SimPar.Re)*(Operators_psi.del4 * Psi.Current)  +...
      + w1; 

  
  
     
 % Boundary Condition 
 N=length(Grid.x)-1;
 S= [1;1-Grid.x(2:N).^2;1];
 RHS(Grid.Walls)=0;  % no velocity on the boundaries ...
 RHS(Grid.TopWall)= -(0.5 * U_lid) ./ S; % except the lid
 
 
 
 % test - remove 
%  RHS =(Operators_psi.del4 * Psi.Current);
end

