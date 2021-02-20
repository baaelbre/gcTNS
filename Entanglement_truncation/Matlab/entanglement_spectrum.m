clear
clc

mu = 1;
nu = 0.4;
D = 2;
% number of iterations: (nmax+1)^D (so can go pretty big as D is small in
% general)
nmax = 30 ;

[Nu, S, Epsilon, xi, entropy, V, alpha] = symplectic_decomposition(D, mu, nu);
V = diag(V);

% brute force way: a nested for loop in which you fill each mode up to a
% certain n_max, and calculate the singular value for each configuration
% every particle in mode i corresponds to a factor xi_i in the Schmidt
% value
vacuum = zeros([D,1]);
states = {vacuum};
schmidts = [prod(1-xi)];
newstate = vacuum;
ctr=2;
for i = 1:D
    for j = 1:length(states)
        newstate = states{j};
        for k = 1:nmax
            states{ctr} = newstate;
            states{ctr}(i) = newstate(i)+1;
            newstate = states{ctr};
            ctr = ctr + 1 ;
        end
    end
end

for i = 2:length(states)
    schmidts(i) = prod(1-xi);
    for j = 1:D
        n = states{i}(j);
        schmidts(i) = schmidts(i)*xi(j)^n;
    end
end

[schmidts, sortIdx] = sort(schmidts, 'descend');
states = states(sortIdx);

schmidt_max = 40;
schmidt_nmbr = linspace(1,schmidt_max,schmidt_max);

ticks = {};
% construct list of state strings
for i = 1:schmidt_max
    ket = states{i};
    str = "";
    str = str + '$$|';
    for j = 1:D
        str = str + num2str(ket(j)) + ', ';
    end
    temp = char(str);
    temp = temp(1:end-2);
    str = string(temp);
    str = str+'>$$';
    ticks{i} = str;
end

set(groot,'defaultAxesTickLabelInterpreter','latex');  
figure();
semilogy(schmidt_nmbr, schmidts(1:schmidt_max),'.', 'MarkerSize', 20);
grid on;
xlabel("state i",'Interpreter','Latex','FontSize',20);
xticks(schmidt_nmbr);
xticklabels(ticks);
ylabel("Schmidt value $$\lambda_i$$",'Interpreter','Latex','FontSize',20);
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$ ', num2str(mu), ', $\nu=$ ', num2str(nu)), 'Interpreter','Latex','FontSize',20);

% matrix elements

T = inv(S');
%a = T(1:D,1:D), %d = T(D+1:2*D,D+1:2*D); both are zero
b = T(1:D,D+1:2*D); c = T(D+1:2*D,1:D);
phi_1 = zeros(schmidt_max, schmidt_max);
% loop over each phi
for M = 1:schmidt_max
    for N = 1:schmidt_max
        R(M,N) = 0;
        Q(M,N) = 0;
        
        % loop over different fields
        for i = 1:D
            phi_i(M,N) = 0;
            phi_i_sq(M,N) = 0;
            pi_i_sq(M,N) = 0;
            
            % equations contain a sum over j
            for j = 1:D
                disp(M), disp(N);
                m = states{M};
                n = states{N};
                phi_ij1 = b(i,j)*n(j)^0.5*eq(m(j), n(j)-1);
                phi_ij2 = b(i,j)*(n(j)+1)^0.5*eq(m(j), n(j)+1);
                
                % the squares: double sums 
                % parts where j = k 
                phisq_ij1 = b(i,j)^2*(n(j)*(n(j)-1))^0.5*eq(m(j), n(j)-2);
                phisq_ij2 = b(i,j)^2*(2*n(j)+1)*eq(m(j), n(j));
                phisq_ij3 = b(i,j)^2*((n(j)+2)*(n(j)+1))^0.5*eq(m(j), n(j)+2);
                pisq_ij1 = c(i,j)^2*(n(j)*(n(j)-1))^0.5*eq(m(j), n(j)-2);
                pisq_ij2 = c(i,j)^2*(2*n(j)+1)*eq(m(j), n(j));
                pisq_ij3 = c(i,j)^2*((n(j)+2)*(n(j)+1))^0.5*eq(m(j), n(j)+2);
                for l = 1:D
                    if l == j
                        phi_ij1 = phi_ij1; % do nothing
                    else
                        phi_ij1 = phi_ij1*eq(m(l),n(l));
                        phi_ij2 = phi_ij2*eq(m(l),n(l));
                        phisq_ij1 = phisq_ij1*eq(m(l),n(l));
                        phisq_ij2 = phisq_ij2*eq(m(l),n(l));
                        phisq_ij3 = phisq_ij3*eq(m(l),n(l));
                        pisq_ij1 = pisq_ij1*eq(m(l),n(l));
                        pisq_ij2 = pisq_ij2*eq(m(l),n(l));
                        phisq_ij3 = phisq_ij3*eq(m(l),n(l));
                    end
                end
                
                phi_i(M,N) = phi_i(M,N) - phi_ij1 + phi_ij2;
                phi_i_sq(M,N) = phi_i_sq(M,N) - phisq_ij1 + phisq_ij2 - phisq_ij3;
                pi_i_sq(M,N) = phi_i_sq(M,N) + pisq_ij1 + pisq_ij2 + pisq_ij3;
                
                % parts where j != k
                for k = 1:D
                    phisq_ij4 = b(i,j)*b(i,k)*(n(j)*n(k))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1);
                    phisq_ij5 = b(i,j)*b(i,k)*((n(j)+1)*n(k))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)-1);
                    phisq_ij6 = b(i,j)*b(i,k)*(n(j)*(n(k)+1))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1);
                    phisq_ij7 = b(i,j)*b(i,k)*((n(j)+1)*n(k))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)+1);
                    pisq_ij4 = c(i,j)*b(i,k)*(n(j)*n(k))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1);
                    pisq_ij5 = c(i,j)*b(i,k)*((n(j)+1)*n(k))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)-1);
                    pisq_ij6 = c(i,j)*b(i,k)*(n(j)*(n(k)+1))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1);
                    pisq_ij7 = c(i,j)*b(i,k)*((n(j)+1)*n(k))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)+1);
                    
                    for l = 1:D
                        if l == j || l == k
                            phi_ij1 = phi_ij1;
                        else
                            phisq_ij4 = phisq_ij4*eq(m(l),n(l));
                            phisq_ij5 = phisq_ij5*eq(m(l),n(l));
                            phisq_ij6 = phisq_ij6*eq(m(l),n(l));
                            phisq_ij7 = phisq_ij7*eq(m(l),n(l));
                            pisq_ij4 = pisq_ij4*eq(m(l),n(l));
                            pisq_ij5 = pisq_ij5*eq(m(l),n(l));
                            pisq_ij6 = pisq_ij6*eq(m(l),n(l));
                            phisq_ij7 = phisq_ij7*eq(m(l),n(l));
                        end
                    end
                    phi_i_sq(M,N) = phi_i_sq(M,N) - phisq_ij4 + phisq_ij5 + phisq_ij6 - phisq_ij7;
                    pi_i_sq(M,N) = phi_i_sq(M,N) + pisq_ij4 + pisq_ij5 + pisq_ij6 + pisq_ij7;
                end
                
                
                
            end
            phi_i(M,N) = alpha(i)/2^0.5*phi_i(M,N); % *1i nog
            phi_i_sq(M,N) = V(i)/2*phi_i_sq(M,N);
            pi_i_sq(M,N) = pi_i_sq(M,N)/2;
        end
        R(M,N) = R(M,N) + phi_i(M,N);
        Q(M,N) = Q(M,N) - phi_i_sq(M,N)/2 - pi_i_sq(M,N)/2;
    end
    
end


%%%% ENERGY OF THE CMPS %%%%
rhoL = FindSSL(Q, R, schmidt_max);
rhoR = FindSSR(Q, R, schmidt_max);
rhoR = rhoR/trace(rhoL*rhoR);
QR = Q*R - R*Q;
R2 = R*R;
ekin = trace(rhoL * QR * rhoR * QR');
density = trace(rhoL * R * rhoR * R');
epot = mu * density;
epair = -nu * trace(rhoL*R2*rhoR + rhoL*rhoR*R2');
    %eint = c * trace(rhoL * R2 * rhoR * R2');

e = ekin + epot + epair;
