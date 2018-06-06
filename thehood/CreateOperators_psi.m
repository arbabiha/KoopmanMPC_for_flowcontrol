function [ Operators_psi ] = CreateOperators_psi( D,I)
%CREATEOPERATORS_q creates the Operator acting on psi(x)
Operators_psi.Dx = kron(D,I);
Operators_psi.Dy = kron(I,D);
D2= D^2; D4=D2^2;
Operators_psi.del2 = kron(D2,I)+kron(I,D2);

Operators_psi.del4 = kron(D4,I)+kron(I,D4)+2*kron(D2,I)*kron(I,D2);



end