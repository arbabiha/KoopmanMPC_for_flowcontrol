function Q = getCostMatrix( N )
%COSTCAVITY computes the objective functionction for cavity flow

% cavity grid and PDE operators
persistent Grid Operators



if isempty(Grid)
    [Grid]=CollocationGrid_q(N);    % the computational grid
    [~,Grid.wc] = clencurt(N);Grid.W = kron(Grid.wc,Grid.wc);   % integration weights
     disp('grid created')
end

if isempty(Operators)
    [ Operators ] = CreateOperators_psi( Grid.D,eye(N+1));
    disp('operators for psi created')
end

Dx = sparse(Operators.Dx);
Dy = sparse(Operators.Dy);
W = Grid.W;
WW = [W(:) ; W(:)];
Wdiag = spdiags(WW,0,numel(WW),numel(WW));
Q = [Dx' Dy']*Wdiag*[Dx;Dy];


return

KE = @(psi)  Grid.W* abs((Operators.Dy * psi).^2) +  Grid.W* abs((Operators.Dx * psi).^2);

% 1.Kinetic Energy of difference with the reference state; - convex
C1 = KE(STATE-Ref_STATE);
C1 = C1*1e4;

% 2. the power consumption on top
du_dy = Operators.Dy*(Operators.Dy * STATE);
du_dy_top = du_dy(Grid.TopWall);
C2 = sum( (INPUT .* Grid.wc(:)) .*du_dy_top);
C2= C2/10;

%3. total drag on the floor of cavity
C3 = Grid.wc * du_dy(Grid.BottomWall);
C3=C3/10;

% 4. kietic energy of top lid
C4 = sum( (INPUT .* Grid.wc(:)) .*INPUT);





% a is the weight matrix of the following candidates

c = [C1,C2,C3,C4]*a(:);



end

