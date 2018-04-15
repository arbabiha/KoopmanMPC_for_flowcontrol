%Forms a block diagonal matrix by N times repeating the matrix A
function X = bdiag(A,N)
m = size(A,1);
n = size(A,2);
X = sparse(m*N, n*N);
for i = 1:N
    X( (i-1)*m+1:i*m, (i-1)*n+1:i*n ) = A;
end