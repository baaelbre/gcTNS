function [states, schmidts] = boson_sampling(D,mu,nu, SV_min)

[Nu, S, Epsilon, xi, entropy, V, alpha, e_gctns] = symplectic_decomposition(D, mu, nu);
% lower bound for the Schmidt values
%SV_min = 0.0000000001;
% normalization factor
K = prod(1-xi);
E_max = log(K)- log(SV_min);
states = {{zeros([D,1])}};
schmidts = [K];
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
            newstate = states{ptle_nmbr}{k};
            
            newstate(i) = newstate(i) + 1;
         
            E = sum(newstate.*Epsilon);
            
            % check if the state is not already in the nptle_sector
            % ismember gives logival array, if all these are 1 it means the
            % state is already in there
            if E <= E_max 
                if isempty(nptle_sector)
                    nptle_sector{1} = newstate;
                    schmidts(end+1) = K*exp(-E);
                else if ~ismember(newstate.', [nptle_sector{:}].', 'rows')
                    nptle_sector{end+1} = newstate;
                    schmidts(end+1) = K*exp(-E);
                    end
                end
            end
        end
    end
    states{end+1} = nptle_sector;
end

% concatenate all particle number sectors
states = [states{:}];
% ordering
[schmidts, sortIdx] = sort(schmidts, 'descend');
states = states(sortIdx);



        

