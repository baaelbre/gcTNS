clear
clc

mu = 1;
nu = 0.4;
D = 3;

% dlmread reads data into a matrix, strcat concatenates a string
data = dlmread(strcat('../../Variational_optimization/Matlab/nonrelboson_real/results/',num2str(mu),'_',num2str(nu),'_', num2str(D), '_1'));
% construct coupling matrix 
V = diag(data(1:D));
alpha = data(D+1:2*D);
M = [[V,-kron(alpha,alpha')];[-kron(alpha',alpha),V]];
Omega = M^(1/2); %2Dx2D matrices
Omega_inv = M^(-1/2);

% RDM: only keep upper left blocks of Omega(_inv) in covariance matrix
X = Omega_inv(1:D, 1:D);
P = Omega(1:D, 1:D);
gamma = blkdiag(X,P);

% actual construction of the symplectic matrix that diagonalizes gamma
XP = X*P;
PX = P*X;
[F, Nu] = eig(XP);
[E, Mu] = eig(PX);
S = [[zeros(D),E*Nu^(1/2)];[F*Nu^(1/2),zeros(D)]];

% quantum information concepts
Nu = diag(Nu);
entropy = sum((Nu+1/2*ones([D,1])).*log(Nu+1/2*ones([D,1])) - (Nu-1/2*ones([D,1])).*log(Nu-1/2*ones([D,1])));
% Schmidt values through Tayloring RDM
nmax = 5
Epsilon = log(1+2*Nu)
xi = exp(-Epsilon)
% xi = e^-eppsilon and schmidt values are (1-xi_1)...(1-xi_D)xi_1^n_1
schmidts = 0;
% xi_2^n_2... xi_D^n_D
Number = linspace(1,D,D);
disp(xi)
figure();

plot(Number, Epsilon, '-x');
xlabel("Mode i $(1->D)$",'Interpreter','Latex','FontSize',20);
ylabel("Entanglement energy $\epsilon_i$",'Interpreter','Latex','FontSize',20);
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$', num2str(mu), ', $\nu=$', num2str(nu)), 'Interpreter','Latex','FontSize',20);
%legend(l,'Interpreter','Latex','FontSize',15,'Location','northeast','NumColumns',2)

plot(Number, xi, '-x');
xlabel("Mode i $(1->D)$",'Interpreter','Latex','FontSize',20);
ylabel(" $e^{-\epsilon_i}$",'Interpreter','Latex','FontSize',20);
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$', num2str(mu), ', $\nu=$', num2str(nu)), 'Interpreter','Latex','FontSize',20);


% matrix elements in countable eigenbasis of RDM (number of eigenstates
% corresponds to bond dimension of corresponding cMPS (with finite Hilbert
% space))
