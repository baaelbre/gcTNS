function [e,grad] = Fobj(U,mu,nu)
%Objective function for GcTNS energy density optimization.

%Evaluation energy density.
D = length(U)/2;
V = diag(U(1:D));
alpha = U(D+1:2*D)';
M = [[V,-kron(alpha,alpha')];[-kron(alpha',alpha),V]];
Mp = M^(1/2);
Mm = M^(-1/2);
X = -1/2*Mp(1:D,D+1:end)+mu/2*Mm(1:D,D+1:end)-nu*Mm(1:D,1:D);
e = real(alpha'*X*alpha);

%Determination of the analytical gradient.
grad = zeros(1,2*D);

%Derivatives of e wrt V.
for i = 1:D
    grad(i) = alpha'*(gradX(M,mu,nu,i,i)+gradX(M,mu,nu,i+D,i+D))*alpha;
end

%Derivatives of e wrt alpha.
grad(D+1:end) = 2*(X*alpha);
for i = 1:D
    ei = zeros(D,1);
    ei(i) = 1;
    dMalphai = -kron(alpha,ei')-kron(ei,alpha');
    dXdalphai = zeros(D);
    for u = 1:D
        for v = D+1:2*D
            dXdalphai = dXdalphai + gradX(M,mu,nu,u,v)*dMalphai(u,v-D);
        end
    end
    %{
    %check
    alpha2 = alpha;
    alpha2(i) = alpha2(i) + 10^(-7);
    M2 = [[V,-kron(alpha2,alpha2')];[-kron(alpha2',alpha2),V]];
    Mp2 = M2^(1/2);
    Mm2 = M2^(-1/2);
    X2 = -1/2*Mp2(1:D,D+1:end)+mu/2*Mm2(1:D,D+1:end)-nu*Mm2(1:D,1:D);
    T = (X2-X)*10^7;
    %}
    grad(i+D) = grad(i+D)+alpha'*dXdalphai*alpha;
end

%{
%check
e = energy(U,mu,nu);
for i=1:2*D
    U2 = U;
    U2(i) = U2(i) + 10^(-7);
    e2 = energy(U2,mu,nu);
    T(i) = (e2-e)*10^7;
end
%}
end

