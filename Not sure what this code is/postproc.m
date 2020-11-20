clear
clc

data = dlmread('results/data');

mu = 1;
nulist = [0.25,0.4,0.45,0.495,0.5];
Dmax = [2,3,3,4,5];

l = [];
for i = 1:length(nulist)
    nu = nulist(i);
    mask = data(:,2) == nu;
    e_gctns = data(mask,4);
    epsfun = @(p)-4*nu^2./(sqrt((p.^2+mu).^2-4*nu^2)+p.^2+mu);
    e_ex = 1/2/pi*integral(epsfun,0,Inf);
    relerr = 1-e_gctns./e_ex;
    Ds = data(mask,3);
    l = [l;strcat("$\frac{\nu}{\mu}=$\,",num2str(nu))];
    semilogy(Ds,relerr,'.','MarkerSize',20)
    hold on
end

xlabel("$D$",'Interpreter','Latex','FontSize',15)
ylabel("Relative energy density error",'Interpreter','Latex','FontSize',15)
legend(l,'Interpreter','Latex','FontSize',15,'Location','eastoutside')

%{
p = 0:0.0001:50;

figure
eps = epsfun(p);
plot(p,eps)

figure
C00 = -1/2*eps./sqrt((p.^2+mu).^2-4*nu^2);
loglog(p,C00)

figure
C01 = nu./sqrt((p.^2+mu).^2-4*nu^2);
loglog(p,C01)

%}