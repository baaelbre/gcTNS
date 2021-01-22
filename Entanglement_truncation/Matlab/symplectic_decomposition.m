clear
clc

mu = 1;
nu = 0.4;
D = 1;

data = dlmread('../../Variational_optimization/Matlab/nonrelboson_real/results/data');
mask = boolean((data(:,1) == mu).*(data(:,2) == nu).*(data(:,3) == 1));
e_gctns = data(mask,4);
data = dlmread(strcat('../../Variational_optimization/Matlab/nonrelboson_real/results/',num2str(mu),'_',num2str(nu),'_1_1'),',');
V = data(1:D);
alpha = data(D+1:2*D);
