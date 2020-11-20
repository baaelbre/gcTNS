function [e,grad] = Fobj(U)
%Objective function for GcTNS energy density optimization.

%Evaluation energy density.
D = (length(U)-1)/4;
V = diag(U(2:D+1))+1i*diag(U(D+2:2*D+1));
alpha = U(2*D+2:3*D+1)'+1i*U(3*D+2:end)';
alpha0 = U(1);
M = [[V,-alpha*alpha'];[-conj(alpha)*transpose(alpha),conj(V)]];
Mp = M^(1/2);
Mm = M^(-1/2);
Mmm = M^(-1);
[W,Sigma] = eig(M,'vector');
Wm = W^(-1);

p = alpha0*(1+(transpose(alpha)*Mmm(1:D,1:D)+alpha'*Mmm(D+1:end,1:D))*alpha);               %expectation value psi
pp = 1/2*transpose(alpha)*Mm(1:D,1:D)*alpha;                                                %expectation value psi psi
pdp = 1/2*alpha'*Mm(D+1:end,1:D)*alpha;                                                     %expectation value psi^dag psi
eint = (2*pdp^2+abs(pp)^2) + 4*abs(p)^2*pdp + (pp*conj(p)^2+conj(pp)*p^2) + abs(p)^4;       %expectation value Lieb-Liniger interaction energy
ekin = -1/2*alpha'*Mp(D+1:end,1:D)*alpha;                                                   %expectation value kinetic energy
e = ekin+eint;


%Determination of the analytical gradient.
grad = zeros(1,length(U));

%Derivative of e wrt alpha0.
grad(1) = 2/alpha0*(4*abs(p)^2*pdp + (pp*conj(p)^2+conj(pp)*p^2) + 2*abs(p)^4);

%Derivatives of e wrt V.
Qp = kron(sqrt(Sigma),ones(1,2*D))+kron(ones(2*D,1),transpose(sqrt(Sigma)));
Qm = -1./(diag(sqrt(Sigma))*Qp*diag(sqrt(Sigma)));
Qp = 1./Qp;
for i = 1:D
    dMdVR = Wm(:,i)*W(i,:)+Wm(:,i+D)*W(i+D,:);                                              %derivatives of M and powers hereof wrt V
    dMdVI = 1i*(Wm(:,i)*W(i,:)-Wm(:,i+D)*W(i+D,:));
    dMpdVR = W*(dMdVR.*Qp)*W^(-1);
    dMpdVI = W*(dMdVI.*Qp)*W^(-1);
    dMmdVR = W*(dMdVR.*Qm)*W^(-1);
    dMmdVI = W*(dMdVI.*Qm)*W^(-1);
    dMmmdVR = -Mmm(:,i)*Mmm(i,:)-Mmm(:,i+D)*Mmm(i+D,:);
    dMmmdVI = -1i*Mmm(:,i)*Mmm(i,:)+1i*Mmm(:,i+D)*Mmm(i+D,:);
    
    DekinR = -1/2*alpha'*dMpdVR(D+1:end,1:D)*alpha;                                         %derivatives of expectation values wrt V
    DekinI = -1/2*alpha'*dMpdVI(D+1:end,1:D)*alpha;
    DpdpR = 1/2*alpha'*dMmdVR(D+1:end,1:D)*alpha;
    DpdpI = 1/2*alpha'*dMmdVI(D+1:end,1:D)*alpha;
    DpR = alpha0*((transpose(alpha)*dMmmdVR(1:D,1:D)+alpha'*dMmmdVR(D+1:end,1:D))*alpha);
    DpI = alpha0*((transpose(alpha)*dMmmdVI(1:D,1:D)+alpha'*dMmmdVI(D+1:end,1:D))*alpha);
    DppR = 1/2*transpose(alpha)*dMmdVR(1:D,1:D)*alpha;
    DppI = 1/2*transpose(alpha)*dMmdVI(1:D,1:D)*alpha;
    DeintR = 4*(pdp+abs(p)^2)*DpdpR+2*real(conj(pp)*DppR)+4*(abs(p)^2+2*pdp)*real(conj(p)*DpR)+2*real(DppR*conj(p)^2+pp*2*conj(p)*conj(DpR));
    DeintI = 4*(pdp+abs(p)^2)*DpdpI+2*real(conj(pp)*DppI)+4*(abs(p)^2+2*pdp)*real(conj(p)*DpI)+2*real(DppI*conj(p)^2+pp*2*conj(p)*conj(DpI));
    
    grad(i+1) = DekinR+DeintR;
    grad(i+D+1) = DekinI+DeintI;
end

%Derivatives of e wrt alpha.
D_ekin = [-real(Mp(D+1:end,1:D)*alpha);-imag(Mp(D+1:end,1:D)*alpha)];                       %derivatives expectation values (contribution outer alpha's)
D_pdp = [real(Mm(D+1:end,1:D)*alpha);imag(Mm(D+1:end,1:D)*alpha)];
D_pp = [Mm(1:D,1:D)*alpha;1i*Mm(1:D,1:D)*alpha];
D_p = 2*alpha0*[Mmm(1:D,1:D)*alpha+real(Mmm(D+1:end,1:D)*alpha);1i*Mmm(1:D,1:D)*alpha+imag(Mmm(D+1:end,1:D)*alpha)];
for i = 1:D
    DekinR = D_ekin(i);
    DekinI = D_ekin(i+D);
    DpdpR = D_pdp(i);
    DpdpI = D_pdp(i+D);
    DpR = D_p(i);
    DpI = D_p(i+D);
    DppR = D_pp(i);
    DppI = D_pp(i+D);
    
    ei = zeros(D,1);                                                                        %derivatives of M and powers hereof wrt alpha
    ei(i) = 1;
    dMdalphaR = [[zeros(D),-ei*alpha'-alpha*ei'];[-ei*transpose(alpha)-conj(alpha)*ei',zeros(D)]];
    dMdalphaI = [[zeros(D),-1i*ei*alpha'+1i*alpha*ei'];[1i*ei*transpose(alpha)-1i*conj(alpha)*ei',zeros(D)]];
    dMpdalphaR = W*((W^(-1)*dMdalphaR*W).*Qp)*W^(-1);
    dMpdalphaI = W*((W^(-1)*dMdalphaI*W).*Qp)*W^(-1);
    dMmdalphaR = W*((W^(-1)*dMdalphaR*W).*Qm)*W^(-1);
    dMmdalphaI = W*((W^(-1)*dMdalphaI*W).*Qm)*W^(-1);
    dMmmdalphaR = -Mmm*dMdalphaR*Mmm;
    dMmmdalphaI = -Mmm*dMdalphaI*Mmm;
    
    DekinR = DekinR -1/2*alpha'*dMpdalphaR(D+1:end,1:D)*alpha;                              %derivatives expectation values (contribution alpha's in M)
    DekinI = DekinI -1/2*alpha'*dMpdalphaI(D+1:end,1:D)*alpha;
    DpdpR = DpdpR +1/2*alpha'*dMmdalphaR(D+1:end,1:D)*alpha;
    DpdpI = DpdpI +1/2*alpha'*dMmdalphaI(D+1:end,1:D)*alpha;
    DpR = DpR +alpha0*((transpose(alpha)*dMmmdalphaR(1:D,1:D)+alpha'*dMmmdalphaR(D+1:end,1:D))*alpha);
    DpI = DpI +alpha0*((transpose(alpha)*dMmmdalphaI(1:D,1:D)+alpha'*dMmmdalphaI(D+1:end,1:D))*alpha);
    DppR = DppR +1/2*transpose(alpha)*dMmdalphaR(1:D,1:D)*alpha;
    DppI = DppI +1/2*transpose(alpha)*dMmdalphaI(1:D,1:D)*alpha;
    DeintR = 4*(pdp+abs(p)^2)*DpdpR+2*real(conj(pp)*DppR)+4*(abs(p)^2+2*pdp)*real(conj(p)*DpR)+2*real(DppR*conj(p)^2+pp*2*conj(p)*conj(DpR));
    DeintI = 4*(pdp+abs(p)^2)*DpdpI+2*real(conj(pp)*DppI)+4*(abs(p)^2+2*pdp)*real(conj(p)*DpI)+2*real(DppI*conj(p)^2+pp*2*conj(p)*conj(DpI));

    grad(i+2*D+1) = DeintR+DekinR;
    grad(i+3*D+1) = DeintI+DekinI;
end

e = real(e);
grad = real(grad);

%{
%check
T = zeros(1,4*D+1);
for i = 1:4*D+1
    U2 = U;
    U2(i) = U2(i) + 10^(-7);
    V2 = diag(U2(2:D+1))+1i*diag(U2(D+2:2*D+1));
    alpha2 = U2(2*D+2:3*D+1)'+1i*U2(3*D+2:end)';
    alpha02 = U2(1);
    M2 = [[V2,-alpha2*alpha2'];[-conj(alpha2)*transpose(alpha2),conj(V2)]];
    Mp2 = M2^(1/2);
    Mm2 = M2^(-1/2);
    Mmm2 = M2^(-1);

    p2 = alpha02*(1+(transpose(alpha2)*Mmm2(1:D,1:D)+alpha2'*Mmm2(D+1:end,1:D))*alpha2);
    pp2 = 1/2*transpose(alpha2)*Mm2(1:D,1:D)*alpha2;
    pdp2 = 1/2*alpha2'*Mm2(D+1:end,1:D)*alpha2;
    eint2 = (2*pdp2^2+abs(pp2)^2) + 4*abs(p2)^2*pdp2 + (pp2*conj(p2)^2+conj(pp2)*p2^2) + abs(p2)^4;
    ekin2 = -1/2*alpha2'*Mp2(D+1:end,1:D)*alpha2;
    e2 = ekin2+eint2;
    T(i) = (e2-e)*10^7;
end
%}
end

