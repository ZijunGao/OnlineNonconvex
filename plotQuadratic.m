% plot : quadratic objection function
%
% data description:
% learning rate; minimal eigenvalue; number of iterations for escape.

%% plot for Gaussian noise
clear all; close all; clc;

%% data
dat = load('quadraticGauDat.mat'); % Gaussian noise
eta = dat.eta;
gamma = dat.gamma;
recordIterationNumberGau = dat.recordIterationNumber;
 
%dat = load('quadraticUniDat.mat'); % Uniform noise
%recordIterationNumberUni = dat.recordIterationNumber;

%dat = load('quadraticBerDat.mat'); % Bernoulli noise
%recordIterationNumberBer = dat.recordIterationNumber;

%dat = load('quadraticSphDat.mat'); % Spherical Gaussian noise
%recordIterationNumberSph = dat.recordIterationNumber;
 
%% plot for Gaussian noise
figure;
hold on;

set(gca,'FontName','Times'),
set(gca,'FontSize',20),
% grid on;
box on;
axis([-inf inf -inf 7*10^5]);
xlabel('$\eta^{-1}$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel('iteration number','FontName','Times','FontSize',20,'interpreter','latex')
title('with extra Gaussian noise','FontName','Times','FontSize',20,'interpreter','latex');

plot((1./eta),(mean(recordIterationNumberGau(1,:,:),3)), 'Color','r','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberGau(2,:,:),3)), 'Color','m','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberGau(3,:,:),3)), 'Color','b','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberGau(4,:,:),3)), 'Color','k','LineWidth', 1.5);

h = legend('$\gamma =-10$', '$\gamma = -5$', '$\gamma = -2.5$', '$\gamma = -1.25$');
set(h,'Interpreter','latex','FontName','Times','FontSize',20,'Location','Northwest') % font type
legend boxon;

set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]);     

saveas(gcf,['./','quadraticGauFig','.pdf']);













%% plot for uniform noise
clear all; close all; clc;

%% data
%dat = load('quadraticGauDat.mat'); % Gaussian noise
%eta = dat.eta;
%b = dat.b;
%recordIterationNumberGau = dat.recordIterationNumber;
 
dat = load('quadraticUniDat.mat'); % Uniform noise
eta = dat.eta;
gamma = dat.gamma;
recordIterationNumberUni = dat.recordIterationNumber;

%dat = load('quadraticBerDat.mat'); % Bernoulli noise
%recordIterationNumberBer = dat.recordIterationNumber;

%dat = load('quadraticSphDat.mat'); % Spherical Gaussian noise
%recordIterationNumberSph = dat.recordIterationNumber;
 
%% plot for uniform noise
figure;
hold on;

set(gca,'FontName','Times'),
set(gca,'FontSize',20),
% grid on;
box on;
axis([-inf inf -inf 7*10^5]);
xlabel('$\eta^{-1}$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel('iteration number','FontName','Times','FontSize',20,'interpreter','latex')
title('with extra uniform noise','FontName','Times','FontSize',20,'interpreter','latex');

plot((1./eta),(mean(recordIterationNumberUni(1,:,:),3)), 'Color','r','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberUni(2,:,:),3)), 'Color','m','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberUni(3,:,:),3)), 'Color','b','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberUni(4,:,:),3)), 'Color','k','LineWidth', 1.5);

h = legend('$\gamma =-10$', '$\gamma = -5$', '$\gamma = -2.5$', '$\gamma = -1.25$');
set(h,'Interpreter','latex','FontName','Times','FontSize',20,'Location','Northwest') % font type
legend boxon;

set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]);     

saveas(gcf,['./','quadraticUniFig','.pdf']);















%% plot for Bernoulli noise
clear all; close all; clc;

%% data
%dat = load('quadraticGauDat.mat'); % Gaussian noise
%recordIterationNumberGau = dat.recordIterationNumber;
 
%dat = load('quadraticUniDat.mat'); % Uniform noise
%recordIterationNumberUni = dat.recordIterationNumber;

dat = load('quadraticBerDat.mat'); % Bernoulli noise
eta = dat.eta;
gamma = dat.gamma;
recordIterationNumberBer = dat.recordIterationNumber;

%dat = load('quadraticSphDat.mat'); % Spherical Gaussian noise
%recordIterationNumberSph = dat.recordIterationNumber;
 
%% plot for Bernoulli noise
figure;
hold on;

set(gca,'FontName','Times'),
set(gca,'FontSize',20),
% grid on;
box on;
axis([-inf inf -inf 7*10^5]);
xlabel('$\eta^{-1}$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel('iteration number','FontName','Times','FontSize',20,'interpreter','latex')
title('with extra Bernoulli noise','FontName','Times','FontSize',20,'interpreter','latex');

plot((1./eta),(mean(recordIterationNumberBer(1,:,:),3)), 'Color','r','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberBer(2,:,:),3)), 'Color','m','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberBer(3,:,:),3)), 'Color','b','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberBer(4,:,:),3)), 'Color','k','LineWidth', 1.5);

h = legend('$\gamma =-10$', '$\gamma = -5$', '$\gamma = -2.5$', '$\gamma = -1.25$');
set(h,'Interpreter','latex','FontName','Times','FontSize',20,'Location','Northwest') % font type
legend boxon;

set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]);     

saveas(gcf,['./','quadraticBerFig','.pdf']);
















%% plot for spherial Gaussian noise
clear all; close all; clc;

%% data
%dat = load('quadraticGauDat.mat'); % Gaussian noise
%recordIterationNumberGau = dat.recordIterationNumber;
 
%dat = load('quadraticUniDat.mat'); % Uniform noise
%recordIterationNumberUni = dat.recordIterationNumber;

%dat = load('quadraticBerDat.mat'); % Bernoulli noise
%recordIterationNumberBer = dat.recordIterationNumber;

dat = load('quadraticSphDat.mat'); % Spherical Gaussian noise
eta = dat.eta;
gamma = dat.gamma;
recordIterationNumberSph = dat.recordIterationNumber;
 
%% log scale plot for spherical Gaussian noise
figure;
hold on;

set(gca,'FontName','Times'),
set(gca,'FontSize',20),
% grid on;
box on;
axis([-inf inf -inf 7*10^5]);
xlabel('$\eta^{-1}$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel('iteration number','FontName','Times','FontSize',20,'interpreter','latex')
title('with extra uniform-spherical noise','FontName','Times','FontSize',20,'interpreter','latex');

plot((1./eta),(mean(recordIterationNumberSph(1,:,:),3)), 'Color','r','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberSph(2,:,:),3)), 'Color','m','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberSph(3,:,:),3)), 'Color','b','LineWidth', 1.5);
plot((1./eta),(mean(recordIterationNumberSph(4,:,:),3)), 'Color','k','LineWidth', 1.5);

h = legend('$\gamma =-10$', '$\gamma = -5$', '$\gamma = -2.5$', '$\gamma = -1.25$');
set(h,'Interpreter','latex','FontName','Times','FontSize',20,'Location','Northwest') % font type
legend boxon;

set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]);     

saveas(gcf,['./','quadraticSphFig','.pdf']);



    
%% bugs and tips:
% 1. x. * y. -> x. * y
% 2. z = x.^2 *¡¡(cos(y)) * -> z = (cos(y)) * x.^2;




