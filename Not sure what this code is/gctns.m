function [] = gctns(D,mu,nu,varargin)
%Optimize GcTNS for non-relativistiv boson in 1d and D bond fields.

%initialize constrained optimization
fun = @(U)Fobj(U,mu,nu);
con = @(U)constraint(U);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy','MaxFunctionEvaluations',25000,'SpecifyObjectiveGradient',true,'StepTolerance',1e-10,'OptimalityTolerance',5e-09,'MaxIterations',10000);

%initial guess
if isempty(varargin)
    flag = true;
    while flag 
        x0 = rand(1,2*D);
        flag = any(con(x0) > 0);
    end
else
end

%actual optimization
x = fmincon(fun,x0,[],[],[],[],[],[],con,options);

%output
files = dir(strcat('results/',num2str(mu),'_',num2str(nu),'_',num2str(D),'_*'));
z = length(files)+1;
dlmwrite(strcat('results/',num2str(mu),'_',num2str(nu),'_',num2str(D),'_',num2str(z)),x,'precision','%18.15f');
fid = fopen('results/data','a');
fprintf(fid,strcat(num2str(mu),'\t',num2str(nu),'\t',num2str(D),'\t',num2str(fun(x),'%18.15f'),'\n'));
fclose(fid);

end

