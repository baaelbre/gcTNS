clear
clc

mu = 1;
nu = 0.4;
D=1;

data = dlmread('../../Variational_optimization/Matlab/nonrelboson_real/results/data');
mask = boolean((data(:,1) == mu).*(data(:,2) == nu).*(data(:,3) == 1));
e_gctns = data(mask,4);
data = dlmread(strcat('../../Variational_optimization/Matlab/nonrelboson_real/results/',num2str(mu),'_',num2str(nu),'_1_1'),',');
V = data(1);
alpha = data(2);
Omega = sqrt(V^2-alpha^4);
epsfun = @(p)-4*nu^2./(sqrt((p.^2+mu).^2-4*nu^2)+p.^2+mu);
e_ex = 1/2/pi*integral(epsfun,0,Inf);

chimax = 15;
e = zeros(1,chimax);
for chi = 1:chimax
    R = alpha/sqrt(2*Omega)*(diag(sqrt(1:chi-1),1)+diag(sqrt(1:chi-1),-1));
    Q = -(V/2/Omega+ Omega/3)*diag(1:chi) + 1/4*(Omega-V/Omega)*(diag(sqrt((2:chi-1).*(1:chi-2)),2)+diag(sqrt((2:chi-1).*(1:chi-2)),-2));

    rhoL = FindSSL(Q, R, chi);
    rhoR = FindSSR(Q, R, chi);

    rhoR = rhoR/trace(rhoL*rhoR);

    QR = Q*R - R*Q;
    R2 = R*R;
    ekin = trace(rhoL * QR * rhoR * QR');
    density = trace(rhoL * R * rhoR * R');
    epot = mu * density;
    epair = -nu * trace(rhoL*R2*rhoR + rhoL*rhoR*R2');
    %eint = c * trace(rhoL * R2 * rhoR * R2');

    e(chi) = ekin + epot + epair;
end

figure
plot(1:chimax,e,'b-d','LineWidth',2)
hold on
plot(0:chimax,e_gctns*ones(1,chimax+1),'r--','LineWidth',2)
plot(0:chimax,e_ex*ones(1,chimax+1),'k--','LineWidth',2)
hold off
legend(["Truncated gcTNS","Full gcTNS","Exact"],'Interpreter','Latex','FontSize',15)
xlabel("$\chi$",'Interpreter','Latex','FontSize',15)
ylabel("$e(\chi)$",'Interpreter','Latex','FontSize',15)
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$', num2str(mu), ', $\nu=$', num2str(nu)), 'Interpreter','Latex','FontSize',20);