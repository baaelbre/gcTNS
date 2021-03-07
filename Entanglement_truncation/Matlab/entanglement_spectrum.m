clear
clc

mu = 1;
nu = 0.4;
D = 1;
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

celldisp(states)

for i = 2:length(states)
    schmidts(i) = prod(1-xi);
    for j = 1:D
        n = states{i}(j);
        schmidts(i) = schmidts(i)*xi(j)^n;
    end
end

[schmidts, sortIdx] = sort(schmidts, 'descend');
states = states(sortIdx);

schmidt_max = 30;
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
T = S.';
A = T(D+1:2*D,1:D); B = T(1:D,D+1:2*D);
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
                %disp(M), disp(N);
                m = states{M};
                n = states{N};
                %phi
                phi_ij1 = B(i,j)*n(j)^0.5*eq(m(j), n(j)-1)/sqrt(2);
                phi_ij2 = B(i,j)*(n(j)+1)^0.5*eq(m(j), n(j)+1)/sqrt(2);
                for l = 1:D
                    if l == j
                        phi_ij1 = phi_ij1; % do nothing
                    else
                        phi_ij1 = phi_ij1*eq(m(l),n(l));
                        phi_ij2 = phi_ij2*eq(m(l),n(l));
                    end
                end
                phi_i(M,N) = phi_i(M,N) - 1i*phi_ij1 + 1i*phi_ij2;
                

                
                for k = 1:D
                    phisq_ij1 = B(i,j)*B(i,k)*(n(j)*n(k))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1)/2;
                    phisq_ij2 = B(i,j)*B(i,k)*((n(j)+1)*n(k))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)-1)/2;
                    phisq_ij3 = B(i,j)*B(i,k)*((n(j)+1)*(n(k)+1))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1)/2;
                    phisq_ij4 = B(i,j)*B(i,k)*((n(j))*(n(k)+1))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)+1)/2;
                    pisq_ij1 = A(i,j)*A(i,k)*(n(j)*n(k))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1)/2;
                    pisq_ij2 = A(i,j)*A(i,k)*((n(j)+1)*n(k))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)-1)/2;
                    pisq_ij3 = A(i,j)*A(i,k)*((n(j)+1)*(n(k)+1))^0.5*eq(m(j), n(j)-1)*eq(m(k),n(k)-1)/2;
                    pisq_ij4 = A(i,j)*A(i,k)*(n(j)*(n(k)+1))^0.5*eq(m(j), n(j)+1)*eq(m(k),n(k)+1)/2;
                    
                    for l = 1:D
                        if l == j || l == k
                            phi_ij1 = phi_ij1;
                        else
                            phisq_ij1 = phisq_ij1*eq(m(l),n(l));
                            phisq_ij2 = phisq_ij2*eq(m(l),n(l));
                            phisq_ij3 = phisq_ij3*eq(m(l),n(l));
                            phisq_ij4 = phisq_ij4*eq(m(l),n(l));
                            pisq_ij1 = pisq_ij1*eq(m(l),n(l));
                            pisq_ij2 = pisq_ij2*eq(m(l),n(l));
                            pisq_ij3 = pisq_ij3*eq(m(l),n(l));
                            phisq_ij4 = phisq_ij4*eq(m(l),n(l));
                        end
                    end
                    phi_i_sq(M,N) = phi_i_sq(M,N) - phisq_ij1 + phisq_ij2 + phisq_ij3 - phisq_ij4;
                    pi_i_sq(M,N) = phi_i_sq(M,N) + pisq_ij1 + pisq_ij2 + pisq_ij3 + pisq_ij4;
                end
                
                
                
            end
            phi_i(M,N) = alpha(i)*phi_i(M,N); 
            phi_i_sq(M,N) = V(i)*phi_i_sq(M,N);
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
