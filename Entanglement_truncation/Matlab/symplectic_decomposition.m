function [Nu, S, Epsilon, xi, entropy, V, alpha] = symplectic_decomposition(D,mu,nu)
% dlmread reads data into a matrix, strcat concatenates a string
data = dlmread(strcat('../../Variational_optimization/Matlab/nonrelboson_real/results/',num2str(mu),'_',num2str(nu),'_', num2str(D), '_1'));

% construct 2Dx2D coupling matrix 
V = diag(data(1:D));
alpha = data(D+1:2*D);
M = [[V,-kron(alpha,alpha')];[-kron(alpha',alpha),V]];
Omega = M^(1/2); 
Omega_inv = M^(-1/2);

% RDM: only keep upper left blocks of Omega(_inv) in covariance matrix
% gamma
X = Omega_inv(1:D, 1:D);
P = Omega(1:D, 1:D);
gamma = blkdiag(X,P);

% actual construction of the symplectic matrix that diagonalizes gamma
XP = (X*P);
PX = (P*X);
[F, Nu] = eig((XP)^1/2);
[E, Mu] = eig((PX)^1/2);
Nu = abs(Nu);

% sort symplectic eigenvalues such that they descend
Nu = diag(Nu);
[Nu, sortIdx] = sort(Nu, 'descend');
E = E(:,sortIdx);
F = F(:,sortIdx);

% according to eq. (3.22) in thesis Quinten, the eigenvectors should be
% orthonormalized according to
E = orthonormalization(E, X);
F = orthonormalization(F, P);
S = [[zeros(D),E*diag(Nu)^(1/2)];[F*diag(Nu)^(1/2),zeros(D)]];

% quantum information concepts

arr = 1/2*ones([D,1]);
entropy = sum((Nu+arr).*log(Nu+arr) - (Nu-arr).*log(Nu-arr));
disp(Nu);

% Schmidt values through Tayloring RDM
nmax = 5;
Epsilon = log(Nu+1/2)-log(Nu-1/2);
xi = exp(-Epsilon);
% xi = e^-eppsilon and schmidt values are (1-xi_1)...(1-xi_D)xi_1^n_1
schmidts = 0;
% xi_2^n_2... xi_D^n_D

% matrix elements in countable eigenbasis of RDM (number of eigenstates
% corresponds to bond dimension of corresponding cMPS (with finite Hilbert
% space))
