function [GX] = gradX(M,mu,nu,i,j)
%Derivative of X wrt M_ij.

D = length(M)/2;
[U,Sigma] = eig(M);
Mtilde = zeros(2*D);
Mtilde(i,j) = 1;
Mtilde(j,i) = 1;
Z = U'*Mtilde*U;
Sigmatilde = diag(diag(Z));
Utilde = zeros(2*D);
for a = 1:2*D
    for b = a+1:2*D
        Utilde(a,b) = 1/(Sigma(b,b)-Sigma(a,a))*Z(a,b);
    end
end
Utilde = Utilde-Utilde';
Gp = U*(Utilde*Sigma^(1/2)-Sigma^(1/2)*Utilde+1/2*diag(diag(Sigmatilde)./(diag(Sigma).^(1/2))))*U';
Gm = U*(Utilde*Sigma^(-1/2)-Sigma^(-1/2)*Utilde-1/2*diag(diag(Sigmatilde)./(diag(Sigma).^(3/2))))*U';

GX = -1/2*Gp(1:D,D+1:end)+mu/2*Gm(1:D,D+1:end)-nu*Gm(1:D,1:D);

%{
%check
Mp = M^(1/2);
Mm = M^(-1/2);
X = -1/2*Mp(1:D,D+1:end)+mu/2*Mm(1:D,D+1:end)-nu*Mm(1:D,1:D);
M2 = M + Mtilde*10^(-7);
Mp2 = M2^(1/2);
Mm2 = M2^(-1/2);
X2 = -1/2*Mp2(1:D,D+1:end)+mu/2*Mm2(1:D,D+1:end)-nu*Mm2(1:D,1:D);
T = (X2-X)*10^7;
%}
end