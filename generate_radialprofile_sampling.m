% Generate the sampling_*.mat files that describe the radial profile of the emitting state
% This script uses Thomson Scattering data, CoRa and GotoCRM to predict spatial profiles of excited states of He I
% We use inverse transform sampling: for each radial profile, calculate the CDF and its inverse and save it in sampling_*.mat

close all
clear all

%% Settings
doPlot = 0;
save_cdf = 1; 

%emitting_state = 11; % 31P = 11;
TSdata_filename = 'TSdata_OES_TS3';
%TScase = 'half';
useOpacity = 0;

%% Load TS data
load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/', TSdata_filename, '.mat']);
data=data(data.Temperature>=0.5,:); % Remove points at too low Te where CRMs do not work
radius = 1e-3*table2array(data(:,6))-5e-3; % Convert to mm and subtract the observed 5mm position shift between OES and TS data

ne = 1e-6*table2array(data(:,7));
Te = table2array(data(:,9));

%% Run CRMs at all TS positions to obtain radial profiles in cartesian coords
c = crm();
c.settings.methodEIE = 'accurate';
c.settings.enableOpacity = useOpacity;
c.settings.resolve_ions = 0;
c.nHe = 5e19; %(check for specific RAID condition!)
c.settings.THe_K = 1e3; 

% Loop over positions
for i = 1:length(ne) 
    c.run(ne(i),Te(i));
    c.runGoto;
    density_profiles_Cora(:, i) = c.densities(1:19); % Store density vector, keep only HeI states
    density_profiles_Goto(:, i) = c.Goto.densities;
    %n_emitting(i) = c.densities(emitting_state);
end

% % Group data in a cell array: each cell contains the radial profile of a state
% for j = 1:19
% % Normalize and store profiles in cell array
% profiles_Cora{j} = density_profiles_Cora(j,:)/max(density_profiles_Cora(j,:));
% profiles_Goto{j} = density_profiles_Goto(j,:)/max(density_profiles_Goto(j,:));
% end

% Call CRM to provide n_emitting
%[radius, n_emitting] = CRM_emissivity_radialprofile(emitting_state,filename,1,TScase,useOpacity); % Radial profile of 31p density based on TS data


%% Calculate CDF^-1 for all profiles and collect them in output struct
% We will have to treat r>0 and r<0 positions separately: 'half' and 'neghalf' profiles 
half_inds = find(radius>=0);
neghalf_inds = find(radius<=0);
cdf_val_save = linspace(0,1,1000); % ~arbitrary: select resolution of output invcdf vector

parfor emitting_state = 1:19
invcdf_Cora_half{emitting_state} = create_invcdf(radius(half_inds), density_profiles_Cora(emitting_state,half_inds),cdf_val_save,doPlot);
invcdf_Cora_neghalf{emitting_state} = create_invcdf(-1*radius(neghalf_inds), density_profiles_Cora(emitting_state,neghalf_inds),cdf_val_save,doPlot);
invcdf_Goto_half{emitting_state} = create_invcdf(radius(half_inds), density_profiles_Goto(emitting_state,half_inds),cdf_val_save,doPlot);
invcdf_Goto_neghalf{emitting_state} = create_invcdf(-1*radius(neghalf_inds), density_profiles_Goto(emitting_state,neghalf_inds),cdf_val_save,doPlot);
end

% Create output struct
RTdata = struct();
RTdata.date = c.date;
RTdata.settings = c.settings;
RTdata.settings.nHe = c.nHe;
RTdata.rateMatrixSE = c.processes.rateMatrixSE;
RTdata.cdf_values = cdf_val_save;
RTdata.invcdf_Cora_half = invcdf_Cora_half;
RTdata.invcdf_Cora_neghalf = invcdf_Cora_neghalf;
RTdata.invcdf_Goto_half = invcdf_Goto_half;
RTdata.invcdf_Goto_neghalf = invcdf_Goto_neghalf;

% Save the struct in a .mat file
OPtext = [];
if useOpacity 
    OPtext = '_OP';
end

if save_cdf
    %cdfval_save = linspace(0,1,1e3);
    %invcdfval_save = invcdf_fun(cdfval_save);
    save(['RTdata_sampling_' TSdata_filename(8:end) OPtext '.mat'],'RTdata');
end


%% Function definition - Calculation of CDF^-1
% Convert profiles to cylindrical coords and calculate their CDF^-1

function invcdf_vect = create_invcdf(radius,density_profile,cdf_values_save,doPlot)
    % Treat inputs: turn into column vectors, needed for interp1
    radius = radius(:);
    density_profile = density_profile(:);

    % Calculate PDF in cylindrical coordinates
    profile = @(r) r.*interp1(radius,density_profile,r); % Multiply by jacobian r
    normalization = integral(profile,0,max(radius));
    pdf = @(r) profile(r)/normalization;

    % CDF(r) is the integral of the PDF from 0 to r
    rlims = [0 max(radius)];
    cdf_fun = @(r) cdf(pdf,r,rlims);



    % Invert the CDF to obtain CDF^-1
    invcdf_fun = @(cdfval) invcdf(cdf_fun,cdfval,max(radius)/2);
    invcdf_vect = invcdf_fun(cdf_values_save);


% For testing: plot results to check that random samples are distributed according to pdf
    if doPlot
        % Plot CDF
        figure; hold on;
        rplot = linspace(-0.01,1.3*max(radius),100);
        plot(rplot,cdf(pdf,rplot,rlims))
        xlabel('Radial position (mm)')
        ylabel('CDF')
        
        % Plot CDF^-1
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
        hist = histogram(rpos,'NumBins',ceil(sqrt(Nrolls)),'Normalization','pdf'); % 'Normalization' Divides values by Nrolls for comparison to PDF
        plot(rplot,pdf(rplot))
        xlabel('Radial position (m)')
        title('Random sampling reproduces the PDF as expected')
        legend('Histogram of radius samples (Normalized)', 'PDF')
    end

% Helper functions to calculate CDF and CDF^-1
    function cdfval = cdf(pdf,r,rlims)
        % For loops are needed for integral() and fzero() to handle multiple values
        % rlims needed to extend to any input, needed for fzero
        for i=1:length(r)
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
        % Goal: invert CDF == find the radius corresponding to a given value of the CDF
        for i = 1:length(cdfval)
            fct = @(r) cdf(r) - cdfval(i)- 1e-6*r_start; % Translate the CDF and find the zero
            invcdfval(i) = fzero(fct, r_start);
        end
    end
end

%% Function definition - CRM radial profiles
% Get emissivity of plasma column from TS data + CRM
% Returns normalized profile

% function [radius, profile] = CRM_emissivity_radialprofile(TSfilename,doPlot,TScase,useOpacity)
% if nargin < 2 doPlot = 0; end
% if nargin < 3 TScase = 'half'; end
% if nargin < 4 useOpacity = 0; end
% 
% % Get Thomson scattering data for ne,Te inputs
% 
% switch TScase    
%     %case 'full'
%     %    load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/f', TSfilename, '.mat']);
%     case 'half'
%         load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/', TSfilename, '.mat']);
%         data=data(data.Radius>0,:); % Remove point at zero since current TS data is shifted by exactly that amount
%         data=data(data.Temperature>=0.5,:);
%         radius = 1e-3*table2array(data(:,6))-5e-3; % Subtract TS shift
%     case 'neghalf'
%         load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/', TSfilename, '.mat']);
%         data=data(data.Radius<=0,:); % Remove point at zero since current TS data is shifted by exactly that amount
%         data=data(data.Temperature>=0.5,:);
%         radius = -1e-3*table2array(data(:,6)); % Subtract TS shift
% end
% ne = 1e-6*table2array(data(:,7));
% Te = table2array(data(:,9));
% %radius = 1e-3*table2array(data(:,6))-5e-3; % Subtract TS shift
% 
% c = crm();
% c.settings.methodEIE = 'accurate';
% c.settings.enableOpacity = useOpacity;
% c.settings.resolve_ions = 0;
% c.nHe = 5e19; %(check for specific RAID condition!)
% c.settings.THe_K = 1e3; 
% 
% %ne(end)=[];Te(end)=[];radius(end)=[]; %Use if full TS data, last point at Te=0.04 is unstable 
% 
% for i = 1:length(ne)
%     c.run(ne(i),Te(i));
%     n_emitting(i) = c.densities(emitting_state);
% end
% 
% profile = n_emitting./max(n_emitting);
% 
% % Plot the profile and compare to ne, Te
% if doPlot
% figure; hold on; grid on;
% plot(radius,n_emitting/max(n_emitting))
% plot(radius,ne/max(ne))
% plot(radius,Te/max(Te))
% legend('n_{31P}','ne','Te')
% title('Normalized density profiles (?)')
% xlabel('Radius (m)')
% ylabel('Density')
% end
% end

