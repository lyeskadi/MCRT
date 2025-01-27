% Calculate the CDF and its inverse for a given radial profile
% Can be used to sample from in haloMC.m
close all
clear all
doPlot = 1;
save_cdf = 0;

filename = 'TSdata_OES_TS2';
TScase = 'half';
useOpacity = 1;

% Call CRM to provide n31p
[radius, n31p] = CRM_emissivity_radialprofile(filename,1,TScase,useOpacity); % Radial profile of 31p density based on TS data

% PDF(r) ~ r*n31p to account for cylindrical geometry
%profile = radius.*n31p';

profile = @(r) r.*interp1(radius,n31p',r);
normalization = integral(profile,0,max(radius));
pdf = @(r) profile(r)/normalization;

% CDF(r) is the integral of the PDF from 0 to r
%cdf = @(r) integral(pdf,0,r);
rlims = [0 max(radius)];
cdf_fun = @(r) cdf(pdf,r,rlims);


figure; hold on;
rplot = linspace(-0.01,1.3*max(radius),100);
%plot(rplot,pdf(rplot))
plot(rplot,cdf(pdf,rplot,rlims))

% Invert the cdf to allow inverse sampling from the pdf
invcdf_fun = @(cdfval) invcdf(cdf_fun,cdfval,max(radius)/2);

%% Plot results to check that random samples are distributed according to pdf
if doPlot 
uniform_samples = 0:0.01:1;
figure;
plot(uniform_samples, invcdf_fun(uniform_samples))
title('Radial position samples for an input [0,1]')
xlabel('Input [0,1]')
ylabel('CDF^{-1}, radius of particle')

% Check the result: randomly drawn samples must reproduce the PDF

Nrolls = 1000;
rpos = invcdf_fun(rand(1,Nrolls));
figure; hold on
hist = histogram(rpos,'Normalization','pdf'); % Divides values by Nrolls for comparison to PDF
plot(rplot,pdf(rplot))
xlabel(radius)
title('Random sampling reproduces the PDF as expected')
legend('Histogram of radius samples (Normalized)', 'PDF')
end

%% Save the inverse cdf
OPtext = [];
if useOpacity 
    OPtext = '_OP';
end

if save_cdf
    cdfval_save = linspace(0,1,1e3);
    invcdfval_save = invcdf_fun(cdfval_save);
    save(['sampling_' filename '_' TScase OPtext '.mat'],'cdfval_save','invcdfval_save');
end
%% Function definitions
% For loops are needed for integral() and fzero() to handle multiple values
% rlims needed to extend to any input, needed for fzero
function cdfval = cdf(pdf,r,rlims)
    for i=1:length(r)
        %if ~isnan(pdf(r(i)))
        if r(i)>rlims(2)
            cdfval(i) = 1;
        elseif r(i) < rlims(1)
            cdfval(i) = 0;
        else
            cdfval(i) = integral(pdf,0,r(i));
        end
    end
end

function invcdfval = invcdf(cdf,cdfval,r_start)
% Goal: invert CDF, find the radius corresponding to a given value of the CDF
for i = 1:length(cdfval)
    fct = @(r) cdf(r) - cdfval(i)- 1e-6*r_start; % Translate the CDF and find the zero
    invcdfval(i) = fzero(fct, r_start);
end
end