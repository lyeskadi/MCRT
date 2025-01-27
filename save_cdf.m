% Save CDF^-1 for inverse sampling
close all
clear all

load('TrhoValues.mat');
k0rho = save_k0rho;
T = save_Tmfp;

% figure;
% plot(k0rho,1-T)

cdf_data = [k0rho;T]';

% Invert: find which k0rho gives T=0.2
uniform_samples = [0:0.0001:1];
plot(uniform_samples,invcft(k0rho,T,uniform_samples))
xlabel('Random numbers from uniform distribution [0,1]')
ylabel('k_0\rho, Inverse of the cumulative distribution function')

cftval = uniform_samples;
invcftval = invcft(k0rho,T,uniform_samples);
%save('sampling.mat','cftval','invcftval')

%% Check for consistency
for i=1:10000
    rho(i) = pickTravelDistance(cftval,invcftval);
end

% discard high ones
rho(rho>50) = [];


figure;hold on
histogram(rho,500,'Normalization','pdf')
plot(rho,-dTdrho(1,rho),'x')
title('Sampling the PDF')
legend('Inverse transform samples','Analytical PDF f(\rho)')

%% Plots

%% Function definitions
function k0rho_val = invcft(k0rho,T,cftval)
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