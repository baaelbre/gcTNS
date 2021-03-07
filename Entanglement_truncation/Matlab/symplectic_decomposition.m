function [Nu, S, Epsilon, xi, entropy, V, alpha] = symplectic_decomposition(D,mu,nu)
% dlmread reads data into a matrix, strcat concatenates a string
data = dlmread(strcat('../../Variational_optimization/Matlab/nonrelboson_real/results/',num2str(mu),'_',num2str(nu),'_', num2str(D), '_1'));

% construct 2Dx2D coupling matrix 
V = diag(data(1:D));
alpha = data(D+1:2*D);
M = [[V,-kron(alpha,alpha')];[-kron(alpha',alpha),V]];
Omega = M^0.5; 
Omega_inv = M^(-0.5);

% RDM: only keep upper left blocks of Omega(_inv) in covariance matrix
% gamma
X = Omega_inv(1:D, 1:D)/2;
P = Omega(1:D, 1:D)/2;
gamma = blkdiag(X,P);

% actual construction of the symplectic matrix that diagonalizes gamma
XP = (X*P);
PX = (P*X);
[F, Nu] = eig((XP)^0.5);
[E, Mu] = eig((PX)^0.5);
Nu = abs(Nu);

% sort symplectic eigenvalues such that they descend
Nu = diag(Nu);
[Nu, sortIdx] = sort(Nu, 'descend');
E = E(:,sortIdx);
F = F(:,sortIdx);

% according to eq. (3.22) in thesis Quinten, the eigenvectors should be
% orthonormalized according to
E = orthonormalization(E, X);
F = -orthonormalization(F, P);
S = [[zeros(D),E*diag(Nu)^(1/2)];[F*diag(Nu)^(1/2),zeros(D)]];

% quantum information concepts

arr = 1/2*ones([D,1]);
entropy = sum((Nu+arr).*log(Nu+arr) - (Nu-arr).*log(Nu-arr));
disp(Nu);

% Schmidt values through Tayloring RDM
nmax = 5;
Epsilon = log(Nu+1/2)-log(Nu-1/2);
xi = exp(-Epsilon);

