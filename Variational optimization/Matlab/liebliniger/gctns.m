function [] = gctns(D,gamma,varargin)
%Optimize GcTNS  with D bond fields for Lieb-Liniger model with gamma = c(=1)/rho.

%initial guess
if isempty(varargin)
    flag = true;
    while flag 
        x0 = rand(1,4*D+1);
        [c,ceq] = constraint_init(x0,gamma);
        flag = any([c;ceq] > 0);
    end
    x0(1) = -ceq;
else
    x0 = varargin;
end

%initialize constrained optimization
fun = @(U)Fobj(U);
con = @(U)constraint(U,gamma);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy','MaxFunctionEvaluations',50000,'SpecifyObjectiveGradient',true,'StepTolerance',1e-10,'OptimalityTolerance',1e-06,'MaxIterations',10000);

%actual optimization
[x,fval,exitflag] = fmincon(fun,x0,[],[],[],[],[],[],con,options);

%output
files = dir(strcat('results/',num2str(D),'_*'));
z = length(files)+1;
dlmwrite(strcat('results/',num2str(D),'_',num2str(z)),x,'precision','%18.15f');
fid = fopen('results/data','a');
fprintf(fid,strcat(num2str(D),'\t',num2str(gamma),'\t',num2str(fval*gamma^3,'%18.15f'),'\t',num2str(exitflag),'\n'));
fclose(fid);

end

