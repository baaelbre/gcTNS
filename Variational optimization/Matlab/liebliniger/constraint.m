function [c,ceq] = constraint(U,gamma)
%Constraint evaluation by determination of eigenvalues and rho.

D = (length(U)-1)/4;
V = diag(U(2:D+1))+1i*diag(U(D+2:2*D+1));
alpha = U(2*D+2:3*D+1)'+1i*U(3*D+2:end)';
alpha0 = U(1);
M = [[V,-alpha*alpha'];[-conj(alpha)*transpose(alpha),conj(V)]];
Mm = M^(-1/2);
Mmm = M^(-1);

c = -real(eig(M));
ceq = real(1/2*alpha'*Mm(D+1:end,1:D)*alpha + abs(alpha0*(1+(transpose(alpha)*Mmm(1:D,1:D)+alpha'*Mmm(D+1:end,1:D))*alpha))^2-1/gamma);

end