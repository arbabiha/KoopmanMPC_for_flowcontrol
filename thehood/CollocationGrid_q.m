function [Grid]=CollocationGrid_q(N)

  [Grid.D,Grid.x] = cheb(N);
  [xx,yy]=meshgrid(Grid.x,Grid.x); Grid.xx=xx(:);Grid.yy=yy(:); 
  
% Index positions for walls
  Grid.Walls = find(abs(xx)==1 | abs(yy)==1);
  Grid.LeftWall = find(xx==-1); Grid.RightWall=find(xx==1);
  Grid.BottomWall=find(yy==-1); Grid.TopWall=find(yy==1);
  
  
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