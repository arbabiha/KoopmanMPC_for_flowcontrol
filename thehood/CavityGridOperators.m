function [ Grid,varargout ] = CavityGridOperators( N )
%CAVITYGRIDOPERATORS generates the grid and operators
addpath('C:\Users\harbabi\Documents\MATLAB\CavityFlow\FlowSimulation\modules')
[Grid]=CollocationGrid_q(N);    % grid coordiuantes and indices

if nargout==2
    varargout{1}= CreateOperators_psi( Grid.D,eye(N+1));  %  derivative matrices
end

[~,Grid.wc] = clencurt(N);
 Grid.W = kron(Grid.wc,Grid.wc);   % integration coefficients for computation of Kinetic Energy
 rmpath('C:\Users\harbabi\Documents\MATLAB\CavityFlow\FlowSimulation\modules')

end

function [Grid]=CollocationGrid_q(N)

  [Grid.D,Grid.x] = cheb(N);
  [xx,yy]=meshgrid(Grid.x,Grid.x); Grid.xx=xx(:);Grid.yy=yy(:); 
  
  
  
% Index positions for walls and interior
  Grid.Walls = find(abs(xx)==1 | abs(yy)==1);
  Grid.LeftWall = find(xx==-1); Grid.RightWall=find(xx==1);
  Grid.BottomWall=find(yy==-1); Grid.TopWall=find(yy==1);
  Grid.Interior = setdiff((1:length(Grid.xx)),Grid.Walls);

  
  % operators
  I = eye(N+1);
  Grid.D2 = (diag(1-Grid.x.^2)*Grid.D^2 - 4*diag(Grid.x)*Grid.D - 2*I);
  Grid.D4 = (diag(1-Grid.x.^2)*Grid.D^4 - 8*diag(Grid.x)*Grid.D^3 - 12*Grid.D^2); % Trefethen style
  

     % Propagation function
   Grid.S0 = (1-Grid.xx.^2).*(1-Grid.yy.^2);
   Grid.S1 = Grid.S0;
   Grid.S1(Grid.Walls)=1;
  
   % plot on cheb grid
   Grid.cplot = @(X) contourf(reshape(xx,N+1,N+1),reshape(yy,N+1,N+1),reshape(X,N+1,N+1),50,'LineStyle','None');

end

  function [x,w] = clencurt(N)
  theta = pi*(0:N)'/N; x = cos(theta);
  w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
  if mod(N,2)==0 
    w(1) = 1/(N^2-1); w(N+1) = w(1);
    for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
  else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
  end
  w(ii) = 2*v/N;
  end
  
 function [ Operators_psi ] = CreateOperators_psi( D,I)
%CREATEOPERATORS_q creates the Operator acting on psi(x)
Operators_psi.Dx = kron(D,I);
Operators_psi.Dy = kron(I,D);
D2= D^2; D4=D2^2;
Operators_psi.del2 = kron(D2,I)+kron(I,D2);

Operators_psi.del4 = kron(D4,I)+kron(I,D4)+2*kron(D2,I)*kron(I,D2);



end