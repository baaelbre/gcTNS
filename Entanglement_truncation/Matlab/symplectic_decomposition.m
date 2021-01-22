mu = 1;
nu = 0.4;
D = 2;

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
% and now X != P^-1 
disp(X), disp(P), disp(gamma);
