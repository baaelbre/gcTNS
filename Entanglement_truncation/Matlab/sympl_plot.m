clear
clc

mu = 1;
nu = 0.4;
D = 3;

[Nu, S, Epsilon, xi, entropy] = symplectic_decomposition(D, mu, nu);
Number = linspace(1,D,D);

figure();
plot(Number, Epsilon, '-x');
xlabel("Mode i $(1->D)$",'Interpreter','Latex','FontSize',20);
ylabel("Entanglement energy $\epsilon_i$",'Interpreter','Latex','FontSize',20);
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$', num2str(mu), ', $\nu=$', num2str(nu)), 'Interpreter','Latex','FontSize',20);


figure();
plot(Number, xi, '-x');
xlabel("Mode i $(1->D)$",'Interpreter','Latex','FontSize',20);
ylabel(" $e^{-\epsilon_i}$",'Interpreter','Latex','FontSize',20);
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$', num2str(mu), ', $\nu=$', num2str(nu)), 'Interpreter','Latex','FontSize',20);

