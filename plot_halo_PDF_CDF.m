% Plot PDF, CDF, CDF^-1 for inverse sampling
close all
clear all

load('TrhoValues.mat');
k0rho = save_k0rho;
T = save_Tmfp;

%s=load('sampling.mat');

% figure;
% plot(k0rho,1-T)

cdf_data = [k0rho;T]';

uniform_samples = [0:0.00005:1];
cdfval = uniform_samples;
invcdfval = invcdf(k0rho,T,uniform_samples);

% Check for consistency
for i=1:length(uniform_samples)
    rho(i) = pickTravelDistance(cdfval,invcdfval);
end

% discard high ones
rho(rho>10) = [];

%% Plots
PDF = @(k0rhovar) -dTdrho(1,k0rhovar);
CDF = @(k0rhovar) 1-interp1(k0rho,T,k0rhovar);
INVCDF = @(rn) invcdf(k0rho,T,rn); 
% 1. PDF
k0rho_range = logspace(-3,1.2);

% Figure 1: PDF and CDF for f(rho)
figure;hold on; box on;
plot(k0rho_range,PDF(k0rho_range),'Linewidth',2);
plot(k0rho_range,CDF(k0rho_range),'Linewidth',2);

% x, y, width, height
height = 400;
set(gcf,'Position',[400 1300 1.3*height height])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',15)
ax = gca;
set(gcf,'color','w')

xlabel('Traveled distance $k_0\rho$','Fontsize',19,'Interpreter','latex')
legend('$f(\rho)$ (PDF)','Corresponding CDF','Fontsize',17,'Interpreter','latex')

% figure;
% plot(k0rho_range,PDF(k0rho_range));
% xlabel('Traveled distance $\rho$','Interpreter','latex')
% ylabel('Probability density of absorption at $\rho$','Interpreter','latex')
% title('PDF')
% 
% figure;
% plot(k0rho_range,CDF(k0rho_range));
% xlabel('Traveled distance $rho$','Interpreter','latex')
% ylabel('Cumulative distribution function','Interpreter','latex')
% title('CDF')

% Figure 2: CDF^-1
figure; box on;
uniform_range_restricted = 0.01:0.001:0.99;
plot(uniform_range_restricted,INVCDF(uniform_range_restricted),'Linewidth',2);
set(gcf,'Position',[400 1300 1.2*height height])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',15)
ax = gca;
set(gcf,'color','w')

xlabel('Random number in [0,1[','Fontsize',19,'Interpreter','latex')
ylabel('Traveled distance $k_0\rho$','Fontsize',19,'Interpreter','latex')
legend('CDF$^{-1}$','Fontsize',17,'Interpreter','latex')

%% Figure 3: Verify with histogram
figure;hold on;box on;
histogram(rho,70,'Normalization','pdf')
plot(sort(rho),-dTdrho(1,sort(rho)),'Linewidth',2,'Color',[0.83 0 0])
%title('Sampling the PDF')
%set(gcf,'Position',[400 1300 1.2*height height])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',15)
ax = gca;
set(gcf,'color','w')

xlabel('Traveled distance $k_0\rho$','Fontsize',19,'Interpreter','latex')
ylabel('Samples in bin (normalized)','Fontsize',19,'Interpreter','latex')

legend('Inverse transform samples of $\rho$','Analytical PDF $f(\rho)$','Fontsize',17,'Interpreter','latex')
% Invert: find which k0rho gives T=0.2

% plot(uniform_samples,invcdf(k0rho,T,uniform_samples))
% xlabel('Random numbers from uniform distribution [0,1]')
% ylabel('k_0\rho, Inverse of the cumulative distribution function')

%% Function definitions
function k0rho_val = invcdf(k0rho,T,cftval)
    Tval = 1-cftval;
    for i=1:length(Tval)
    [val, ind] = min(abs(T-Tval(i)));
    k0rho_val(i) = k0rho(ind);
    end
end

function k0rho = pickTravelDistance(cftval,invcftval)
    uniform_sample = rand();
    k0rho = interp1(cftval,invcftval,uniform_sample);
end

function  res = dTdrho(k0,rho)
    for i = 1:length(rho)
        integrand = @(x) exp(-2*x.^2-k0*rho(i)*exp(-x.^2));
        res(i) = -(k0/sqrt(pi))*integral(integrand,-inf,inf);
    end
end