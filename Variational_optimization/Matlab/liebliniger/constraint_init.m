function [c,ceq] = constraint_init(U,gamma)
%Check positivity eigenvalues M and rescale alpha0 (when possible) to have rho = 1/gamma.

D = (length(U)-1)/4;
V = diag(U(2:D+1))+1i*diag(U(D+2:2*D+1));
alpha = U(2*D+2:3*D+1)'+1i*U(3*D+2:end)';
alpha0 = U(1);
M = [[V,-alpha*alpha'];[-conj(alpha)*transpose(alpha),conj(V)]];
Mm = M^(-1/2);
Mmm = M^(-1);

c = -real(eig(M));
ceq = real(1/2*alpha'*Mm(D+1:end,1:D)*alpha + abs(alpha0*(1+(transpose(alpha)*Mmm(1:D,1:D)+alpha'*Mmm(D+1:end,1:D))*alpha))^2-1/gamma);

if abs(ceq) > 10^(-10)
    if 1/gamma-real(1/2*alpha'*Mm(D+1:end,1:D)*alpha) > 0
        ceq = -sqrt(1/gamma-real(1/2*alpha'*Mm(D+1:end,1:D)*alpha))/abs(1+(transpose(alpha)*Mmm(1:D,1:D)+alpha'*Mmm(D+1:end,1:D))*alpha);
    else
        ceq = abs(ceq);
    end
end
end