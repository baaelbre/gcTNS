mu = 1;
nu = 0.4;
D = 3;

[Nu, S, Epsilon, xi, entropy, V, alpha, e_gctns] = symplectic_decomposition(D, mu, nu);
% lower bound for the Schmidt values
SV_min = 0.0000000001;
% normalization factor
K = prod(1-xi);
E_max = log(K)- log(SV_min);
states = {{{zeros([D,1]);K}}};
% obtain the n particle sector with states{n+1}
% max number of particles that we will redistribute
ptle_max = 5;

% loop over all particle sectors
for ptle_nmbr = 1:ptle_max
    nptle_sector = {};
    % each state in the n particle sector, take it and add one extra
    % particle in each mode
    for k = 1: length(states{ptle_nmbr})
        for i= 1:D
            newstate = states{ptle_nmbr}{k}{1};
            
            newstate(i) = newstate(i) + 1;
            disp(newstate);
            E = states{ptle_nmbr}{k}{2} + Epsilon(i);
            disp(E), disp(E_max);
            if E <= E_max
                nptle_sector{end+1} = {newstate; E};
                
            end
        end
    end
    states{end+1} = nptle_sector;
end

% concatenate all particle number sectors
states = [states{:}];
% ordering
schmidts = [];
states_ = {};
for k = 1: length(states)
    schmidts(end+1) = states{k}{2};
    states_{end+1} = states{k}{1};
end
[schmidts, sortIdx] = sort(schmidts);
states = states_(sortIdx);



        

