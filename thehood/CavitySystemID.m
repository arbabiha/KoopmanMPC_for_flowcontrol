function CavitySystemID(DataFileName)

% building the Koopman linear system for cavity flow
% from the data in the file named "DataFileName"
disp('loading EDMD data ...')

% reminder
if ~strcmp(DataFileName,'Cavity_data_4EDMD_0.mat')
        warning(['if yu are using your own data', ...
        'make sure to remove the snapshots that contain NaN or Inf']);
end

load(DataFileName)



% the results are better if we center the data around point close to action
% in the state space
SUBTRACT_MEAN = 1;
if(SUBTRACT_MEAN)
    load('CavityStateLibrary.mat','FixedPoint_Re10k')
    x_mean = FixedPoint_Re10k;
    X = bsxfun(@minus,X,x_mean);
    Y = bsxfun(@minus,Y,x_mean);
end



% the following sensor locations generate the plots in the paper
load('CavitySensorLocations.mat','pindex')


% construct models for various number (k) of spatial measurements
for k = [2,5,10,40,50,100]
    tic
    
    Grid= CavityGridOperators( SimPar.N );
    Grid.Interior_index = setdiff(1:length(Grid.xx),Grid.Walls);
    

    sampind = pindex(1:k);
    % sampind is the index of the measurement point on the Chebyshev grid
    % feel free to replace sampind with some random numbers between 1 and
    % 47^2, 
    
    
    CollectOutput = @(x) x(sampind,:);    
    
 
    
    % location in the cavity domain
    figure(311)
    plot(CollectOutput(Grid.xx),CollectOutput(Grid.yy),'x')
    axis square
    axis([-1,1,-1,1])
    title(['location of sensors in cavity, k=',num2str(k)])
    drawnow
    XX = CollectOutput(X);
    YY = CollectOutput(Y);
    
    
    
    
    % embed the data from each trajectory
    nd = 5; % delay order
    tic
    Xtilde =[]; Ytilde=[]; Utilde=[];

    SimLength=unique(SimLengthVar);
    if length(SimLength)>1
        error('trajectories of different lengtnhs .. adopt')
    end
    for i=1:length(SimLengthVar)

        % consider adopting to different trajectory lengths

        TrajIndex = (i-1)*SimLength + (1:SimLength);    % indexes of each trajectory

        % embed the data
        % this gives upside down - think if ordering matters
        Xe = DelayEmbed(XX(:,TrajIndex) , nd    );
        Ye = DelayEmbed(YY(:,TrajIndex) , nd    );
        Ue = DelayEmbed( U(:,TrajIndex) , nd - 1);


        % add the delay-embedded input as an observable
        Xtilde = [Xtilde,[Xe;Ue(:,1:end-1)]];
        Ytilde = [Ytilde,[Ye;Ue(:,2:end  )]];

        UIndex = (i-1)*SimLength + (nd:SimLength);    % indexes of each trajectory
        Utilde = [Utilde, U(:,UIndex) ];

    end


    % add l-2 norm of recent output and constant as extra observables
    KE_embed = @(y) sum(abs(y((nd-1)*k+1:nd*k,:)).^2,1);
    
    % now lift the data
    X_lift = [Xtilde;KE_embed(Xtilde);ones(1,size(Xtilde,2))];
    Y_lift = [Ytilde;KE_embed(Ytilde);ones(1,size(Ytilde,2))];
    U_lift = Utilde;

    % compute the Koopman linear system 
    % Via EDMD in the normal form
    W = [X_lift;U_lift]*[X_lift;U_lift]';
    V = Y_lift*[X_lift;U_lift]';
    M = V*pinv(W);



    Nlift = size(X_lift,1);
    A = M(:,1:Nlift);
    B = M(:,Nlift+1:end);
    C = zeros(k,size(A,1));
    C(:,(nd-1)*k+1:nd*k)=eye(k);

    
    
    % define a lift function to be used in control:
    % first a simple dely-embed that lines up the last nd columns of given
    % data matrix
    Delay_Embed = @(X,nd) reshape(X(:,end-nd+1:end),nd*size(X,1),1);
    
    % then a function that forms zeta given the state and input matrices
    CreateZeta = @(X,U) [Delay_Embed(CollectOutput(X(:,end-nd+1:end)),nd)...
                            ;Delay_Embed(U(:,1:end-1),nd-1)];
    
    % then the second lifting
    % the last function gets the actual state and input and returns the
    % Koopman linear state
    LiftFun = @(X,U) [CreateZeta(X,U); sum(CollectOutput(X(:,end)).^2);1];
    
    
    
    
    
    modelfilename = ['KoopmanLinSys_Re13_k',num2str(k),'_random'];
    save(modelfilename,'A','B','C','nd','CollectOutput','KE_embed','x_mean','LiftFun','CreateZeta','Delay_Embed')

    fprintf('Koopman model of cavity created using %d measurements \n',k)
    toc
end
