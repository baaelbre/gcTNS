clear
clc

mu = 1;
nu = 0.4;
D = 3;
% number of iterations: (nmax+1)^D (so can go pretty big as D is small in
% general)
nmax = 4 ;

[Nu, S, Epsilon, xi, entropy] = symplectic_decomposition(D, mu, nu);

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
celldisp(states);

schmidt_max = 5;
schmidt_nmbr = linspace(1,schmidt_max,schmidt_max);

ticks = {};
% construct list of state strings
for i = 1:schmidt_max
    ket = states{i};
    str = "";
    str = str + '$$|'
    for j = 1:D
        str = str + num2str(ket(j)) + ', ';
    end
    temp = char(str);
    temp = temp(1:end-2);
    str = string(temp);
    str = str+'>$$';
    ticks{i} = str;
end
celldisp(ticks)

set(groot,'defaultAxesTickLabelInterpreter','latex');  
figure();
plot(schmidt_nmbr, schmidts(1:schmidt_max),'.', 'MarkerSize', 20);
grid on;
xlabel("state i",'Interpreter','Latex','FontSize',20);
xticks(schmidt_nmbr);
xticklabels(ticks);
ylabel("Schmidt value $$\lambda_i$$",'Interpreter','Latex','FontSize',20);
title(strcat('$D=$', num2str(D), ' simulation with', ' $\mu=$ ', num2str(mu), ', $\nu=$ ', num2str(nu)), 'Interpreter','Latex','FontSize',20);

% now calculate for each state the corresponding schmidt value