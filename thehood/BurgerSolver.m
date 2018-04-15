function [X]=BurgerSolver(InitialSTATE,u,SimPar)
N = SimPar.N; 

persistent D D2 

if isempty(D)
    
    D2 = fd(N,3,[0:2*pi/N:(N-1)*2*pi/N],2);
    C(1:N-1,2:N)=eye(N-1); C(N,1)=1;
    D2(1,:)= D2(2,:)*C';
    D2(N,:)= D2(N-1,:)*C;

    D = fd(N,3,[0:2*pi/N:(N-1)*2*pi/N],1);
    D(1,:)= D(2,:)*C';
    D(N,:)= D(N-1,:)*C;
end


% iterate and solve
X = zeros(N,SimPar.N_t);   
X(:,1) = InitialSTATE; % initial condition

for t = 2:SimPar.N_t
    v = X(:,t-1);
    for n = 1
        R = SimPar.dt*(D*v).*v + (eye(N,N) - SimPar.dt*SimPar.nu*D2)*v - X(:,t-1) - SimPar.dt*u; % Residual
        J = 2*SimPar.dt*D*diag(v) + (eye(N,N)-SimPar.dt*SimPar.nu*D2); % Jacobian
        v = v - J\R;
        X(:,t) = v;
    end
end

X = X(:,end);
end

function Dm = fd(N,dd,x,nd)
% finite diff matrix

for n=1:(dd-1)/2
    w=weights(x(n),x(1:dd),nd);
    Dm(n,1:dd)=w(nd+1,:);
end

for n=((dd-1)/2)+1:N-(dd-1)/2
    w=weights(x(n),x(n-((dd-1)/2):n+((dd-1)/2)),nd);
    Dm(n,n-((dd-1)/2):n+((dd-1)/2))=w(nd+1,:);
end
for n=(N-((dd-1)/2))+1:N
    w=weights(x(n),x(N-dd+1:N),nd);
    Dm(n,N-dd+1:N)=w(nd+1,:);
end
end


function w = weights(z,x,m)
% finite diff weights
n=length(x); w=zeros(m+1,n); c1=1; c4=x(1)-z; w(1,1)=1;
for i=2:n
   mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
   for j=1:i-1
      c3=x(i)-x(j);  c2=c2*c3;
      if j==i-1
         w(2:mn,i)=c1*((1:mn-1)'.*w(1:mn-1,i-1)-c5*w(2:mn,i-1))/c2;
         w(1,i)=-c1*c5*w(1,i-1)/c2;
      end
      w(2:mn,j)=(c4*w(2:mn,j)-(1:mn-1)'.*w(1:mn-1,j))/c3;
      w(1,j)=c4*w(1,j)/c3;
   end
   c1=c2;
end
end
