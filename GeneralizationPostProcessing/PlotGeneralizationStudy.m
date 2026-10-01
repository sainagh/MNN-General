clear
clc
close all

%%

data = load('GeneralizationResults.txt');

figure
hold on
plot(data([1 2 4 5 7 8],1),data([1 2 4 5 7 8],2),'.r','MarkerSize',9)
plot([data(1,1) data(4,1) data(7,1)],[mean(data(1:2,2)) mean(data(4:5,2)) mean(data(7:8,2))],'*r','MarkerSize',15,'LineWidth',2.5)
axis([8.5 11.5 0.001 0.1])
box on
set(gca,'LineWidth',2.5)
set(gca,'Xtick',[9 10 11])
yscale('log')
set(gca,'FontSize',18)


figure
hold on
plot(data([1 2 4 5 7 8],1),data([1 2 4 5 7 8],3),'.r','MarkerSize',9)
plot([data(1,1) data(4,1) data(7,1)],[mean(data(1:2,3)) mean(data(4:5,3)) mean(data(7:8,3))],'*r','MarkerSize',15,'LineWidth',2.5)
axis([8.5 11.5 0.001 0.1])
box on
set(gca,'LineWidth',2.5)
set(gca,'Xtick',[9 10 11])
yscale('log')
set(gca,'FontSize',18)



figure
hold on
plot(data([1 2 4 5 7 8],1),data([1 2 4 5 7 8],3)-data([1 2 4 5 7 8],2),'.r','MarkerSize',9)
plot([data(1,1) data(4,1) data(7,1)],[mean(data(1:2,3)-data(1:2,2)) mean(data(4:5,3)-data(4:5,2)) mean(data(7:8,3)-data(7:8,2))],'*r','MarkerSize',15,'LineWidth',2.5)
axis([8.5 11.5 -1e-3 1e-3])
box on
set(gca,'LineWidth',2.5)
set(gca,'Xtick',[9 10 11])
set(gca,'FontSize',18)