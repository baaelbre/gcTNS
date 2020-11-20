clear
clc

data = dlmread('results/data');

Ds = [1,2];
l =[];

for i = 1:length(Ds)
    mask = boolean((data(:,1)==Ds(i)).*(data(:,4)==1));
    dataD = data(mask,:);
    plot(sort(dataD(:,2)),sort(dataD(:,3)),'d-','LineWidth',2)
    hold on
    l = [l,strcat("D = ",num2str(Ds(i)))];
end

xlabel("$\frac{c}{\rho}$",'Interpreter','Latex','FontSize',20)
ylabel("$\frac{e}{\rho^3}$",'Interpreter','Latex','FontSize',20,'Rotation',0)
xlim([0,10])
legend(l,'Interpreter','Latex','FontSize',15,'Location','southeast')
