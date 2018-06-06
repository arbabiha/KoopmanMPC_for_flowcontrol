function H = DelayEmbed(x,nd)
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