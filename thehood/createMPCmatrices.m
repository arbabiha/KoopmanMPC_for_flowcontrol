% Lifted matrices for MPC on horizon N
function [Ab Bb] = createMPCmatrices(A,B,N,includex0)
    if(~exist('includex0','var'))
        includex0 = 0; %first block rows of Ab and Bb corresponding to x0 removed
    end
    
    n = size(A,1);
    Ab = zeros((N+1)*n, n);
    Ab(1:n,:) = eye(n,n);
    for i = 2:N+1
        if( size(A,1) > 1000 )
            i
        end
        %Ab((i-1)*n+1:i*n,:) = A^(i-1);
        Ab((i-1)*n+1:i*n,:) = Ab((i-2)*n+1:(i-1)*n,:)*A;
    end
    if(includex0 == 0)
        Ab = Ab(n+1:end,:); % Seems to take a lot of time if A is huge, possible to speedup
    end
    
    for q = 1:length(B)
        m = size(B{q},2);
        Bb{q} = zeros((N+1)*n, N*m);
        for i = 2:N+1
            Bb{q}((i-1)*n+1:i*n,:) = A * Bb{q}((i-2)*n+1:(i-1)*n,:);
            Bb{q}((i-1)*n+1:n*i,(i-2)*m+1:m*(i-1)) = B{q};
        end
        if(includex0 == 0)
            Bb{q} = Bb{q}(n+1:end,:);
        end
    end


end