function [A,B,C,BuildKoopmanState]= SystemID_via_EDMD(X,Y,U,Ntraj,SimLength,CollectMeasurements,nd)

% pick the data
XX = CollectMeasurements(X);
YY = CollectMeasurements(Y);
nm = size(XX,1); % number of measurements

% embed the data from each trajectory
tic
Xtilde =[]; Ytilde=[]; Utilde=[];

if nd>1

    disp('embedding the data ..')
    for i=1:Ntraj

        % consider adopting to different trajectory lengths

        TrajIndex = (i-1)*SimLength + (1:SimLength);    % indexes of each trajectory

        % embed the data
        Xe = DelayEmbedBlock(XX(:,TrajIndex) , nd    );
        Ye = DelayEmbedBlock(YY(:,TrajIndex) , nd    );
        Ue = DelayEmbedBlock( U(:,TrajIndex) , nd - 1);


        % embed 
        Xtilde = [Xtilde,[Xe;Ue(:,1:end-1)]];
        Ytilde = [Ytilde,[Ye;Ue(:,2:end  )]];

        UIndex = (i-1)*SimLength + (nd:SimLength);    % indexes of each trajectory
        Utilde = [Utilde, U(:,UIndex) ];
        
    end
else
    % in case of no embedding (nd=1)
    Xtilde = XX;
    Ytilde = YY;
    Utilde = U;
    
    
    
end
    
toc
%% Extended Dynamic Mode Decomposition

% l2-norm of the current measurement
KE = @(y) sum(abs(y((nd-1)*nm+1:nd*nm,:)).^2)/nm;

% add l2-norm and constant as extra observables
Xp = [Xtilde;KE(Xtilde);ones(1,size(Xtilde,2))];
Yp = [Ytilde;KE(Ytilde);ones(1,size(Ytilde,2))];
Up = Utilde;
% compute the Koopman linear system 
tic
disp('Running EDMD ...')
W = [Xp;Up]*[Xp;Up]';
V = Yp*[Xp;Up]';
M = V*pinv(W);

Nlift = size(Xp,1);
n = size(XX,1)*nd;
A = M(:,1:Nlift);
B = M(:,Nlift+1:end);

nm= size(XX,1);
C = zeros(nm,size(A,1));
C(:,(nd-1)*nm+1:nd*nm)=eye(nm);
toc
% create a handle for computing the Koopman linear state
Delay_Embed = @(X,nd) reshape(X(:,end-nd+1:end),nd*size(X,1),1);
CreateZeta = @(X,U) [Delay_Embed(CollectMeasurements(X),nd); Delay_Embed(U(:,end-nd+1:end-1),nd-1)];
KE = @(x) sum(abs(x).^2)/length(x);

BuildKoopmanState = @(X,U) [CreateZeta(X,U);KE(CollectMeasurements(X(:,end)));1];




save('KoopmanLinSys_Sparse','A','B','C','nd','BuildKoopmanState','CollectMeasurements','CreateZeta')
end

% delay emebedding of trajectory data
function H = DelayEmbedBlock(x,nd)
% takes a matrix y [size: n * m] and returns its embeded version
% which is H [size: nd*n  * (m-nd+1)]

[n,m]=size(x);
H =zeros(nd*n,m-nd-1);

for j=1:m-nd+1
   % segment of data to be cut out
   Index1 = (j-1)*n + 1;
   Index2 = (j-1)*n + nd*n;
   H(:,j) = x(Index1:Index2);
    
end
end
