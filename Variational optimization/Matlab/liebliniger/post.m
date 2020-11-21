clear
clc

data = dlmread('results/data');

Ds = [1,2];
l = [];

figure('Position',[0,0,800,600])
for i = 1:length(Ds)
    mask = boolean((data(:,2)==Ds(i)).*(data(:,4)==1));
    dataD = data(mask,:);
    plot(sort(dataD(:,1)),sort(dataD(:,3)),'d-','LineWidth',2)
    hold on
    l = [l,strcat("D = ",num2str(Ds(i)))];
end

set(gca,'TickLabelInterpreter','Latex')
set(gca,'FontSize',15)
xl = xlabel("$\frac{c}{\rho}$",'Interpreter','Latex','FontSize',25);
xl.Position(2) = xl.Position(2) +0.2;
yl = ylabel("$\frac{e}{\rho^3}$",'Interpreter','Latex','FontSize',25,'Rotation',0);
yl.Position(1) = yl.Position(1) -0.7;
xlim([0,10])
legend(l,'Interpreter','Latex','FontSize',15,'Location','southeast')
