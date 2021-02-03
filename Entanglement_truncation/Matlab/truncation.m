clear
clc

mu = 1;
nu = 0.5;
D = 3;

[Nu, S, Epsilon, xi, entropy, V, alpha] = symplectic_decomposition(D, mu, nu);
T = inv(S');
a = T(1:D,1:D), b = T(1:D,D+1:2*D), c = T(D+1:2*D,1:D), d = T(D+1:2*D,D+1:2*D);
R = alpha;