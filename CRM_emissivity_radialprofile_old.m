% Get emissivity of plasma column from TS data + CRM
% Used as initial condition for haloMC.m 
%close all
%clear all

% Returns normalized profile
function [radius, profile] = CRM_emissivity_radialprofile_old(emitting_state,TSfilename,doPlot,TScase,useOpacity)
if nargin < 2 doPlot = 0; end
if nargin < 3 TScase = 'half'; end
if nargin < 4 useOpacity = 0; end

% Get Thomson scattering data for ne,Te inputs

switch TScase    
    case 'full'
        load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/f', TSfilename, '.mat']);
    case 'half'
        load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/', TSfilename, '.mat']);
        data=data(data.Radius>0,:); % Remove point at zero since current TS data is shifted by exactly that amount
        data=data(data.Temperature>=0.5,:);
        radius = 1e-3*table2array(data(:,6))-5e-3; % Subtract TS shift
    case 'neghalf'
        load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/', TSfilename, '.mat']);
        data=data(data.Radius<=0,:); % Remove point at zero since current TS data is shifted by exactly that amount
        data=data(data.Temperature>=0.5,:);
        radius = -1e-3*table2array(data(:,6)); % Subtract TS shift
end
ne = 1e-6*table2array(data(:,7));
Te = table2array(data(:,9));
%radius = 1e-3*table2array(data(:,6))-5e-3; % Subtract TS shift

c = crm();
c.settings.methodEIE = 'accurate';
c.settings.enableOpacity = useOpacity;
c.settings.resolve_ions = 0;
c.nHe = 5e19; %(check for specific RAID condition!)
c.settings.THe_K = 1e3; 

%ne(end)=[];Te(end)=[];radius(end)=[]; %Use if full TS data, last point at Te=0.04 is unstable 

for i = 1:length(ne)
    c.run(ne(i),Te(i));
    n_emitting(i) = c.densities(emitting_state);
end

profile = n_emitting./max(n_emitting);

% Plot the profile and compare to ne, Te
if doPlot
figure; hold on; grid on;
plot(radius,n_emitting/max(n_emitting))
plot(radius,ne/max(ne))
plot(radius,Te/max(Te))
legend('n_{31P}','ne','Te')
title('Normalized density profiles (?)')
xlabel('Radius (m)')
ylabel('Density')
end
end

