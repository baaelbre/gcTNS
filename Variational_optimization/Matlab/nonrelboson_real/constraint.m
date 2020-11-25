function [c,ceq] = constraint(U)
%Constraint evaluation by determination of eigenvalues.

D = length(U)/2;
V = diag(U(1:D));
alpha = U(D+1:2*D)';
M = [[V,-kron(alpha,alpha')];[-kron(alpha',alpha),V]];
c = -eig(M);
ceq = [];

end

