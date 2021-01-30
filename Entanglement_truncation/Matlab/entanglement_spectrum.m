clear
clc

mu = 1;
nu = 0.4;
D = 2;
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

schmidt_nmbr = linspace(1,5,5);

% now calculate for each state the corresponding schmidt value