% Run MCRT.m and plot results
close all
clear all
%% General params
L_RAID = 1.5;
Afl = 13372000;
Avuv = 566340000;
branching_ratio = Afl/(Afl+Avuv);

%% Run MCRT
TSfile = 'TS2';
%TScase = 'half';
useOpacity = 0;
Nparticles = 1e5;
nground = 5e19;

% Run MC for each half separately (like OES analysis)
[emission_radius, escapeproba] = MCRT(nground,branching_ratio,TSfile,Nparticles,'half',useOpacity);
[emission_radius_neg, escapeproba_neg] = MCRT(nground,branching_ratio,TSfile,Nparticles,'neghalf',useOpacity);

emission_radius(isnan(emission_radius)) = [];
emission_radius_neg(isnan(emission_radius_neg)) = [];
%% Tally points
filename = ['TSdata_OES_' TSfile];
N_initialemission = 1e6; %floor(branching_ratio/(1-branching_ratio)*Nparticles); % Seems to be inaccurate (was due to low number of samples)
points_initial_emission= haloMC_initialize_plasmacolumn(filename,N_initialemission,L_RAID,'half',useOpacity);%initial_emission_positions; %generateCylinderPoints(R_EXC,L_RAID,N_initialemission);
points_initial_emission_neg= haloMC_initialize_plasmacolumn(filename,N_initialemission,L_RAID,'neghalf',useOpacity);

radius_iem = sqrt(points_initial_emission(1,:).^2+points_initial_emission(2,:).^2);
radius_iem_neg = sqrt(points_initial_emission_neg(1,:).^2+points_initial_emission_neg(2,:).^2);
rplot=1:0.01:100;

radialprofile_halo = radial_profile(emission_radius,rplot,80); % Nbins
radialprofile_halo = radialprofile_halo/max(radialprofile_halo);
radialprofile_halo_neg = radial_profile(emission_radius_neg,rplot,80); % Nbins
radialprofile_halo_neg = radialprofile_halo_neg/max(radialprofile_halo_neg);

radialprofile_init = radial_profile(radius_iem,rplot,100);
radialprofile_init = radialprofile_init/max(radialprofile_init);
radialprofile_init_neg = radial_profile(radius_iem_neg,rplot,100);
radialprofile_init_neg = radialprofile_init_neg/max(radialprofile_init_neg);
%escapeproba = length(wall_absorption_pos)/Nparticles;
%radialprofile_RT = branching_ratio*radialprofile_init + (1-branching_ratio)*radialprofile_halo; % escapeproba %@(r) branching_ratio*radialprofile_init(r) + (1-branching_ratio-escapeproba)*radialprofile_halo(r);
%radialprofile_RT = branching_ratio*radialprofile_init + (1-branching_ratio-escapeproba)*radialprofile_halo;
%radialprofile_RT = radialprofile_halo; %radialprofile_RT/max(radialprofile_RT);

%if TScase == "neghalf"
    rplot_neg = -rplot;
%end

rplot = [flip(rplot_neg) rplot];
radialprofile_RT = [flip(radialprofile_halo_neg) radialprofile_halo];
radialprofile_init = [flip(radialprofile_init_neg) radialprofile_init];

%% Plot radial profile of emission line

% Define colors
% palatinateblue =[0.15, 0.23, 0.89];
% persianblue=[0.11, 0.22, 0.73];
% rosso = [0.83 0 0];
% Azure: 0 0.5 1
% Royal azure [0 0.22 0.66]
% hfacecolor, lightred = [0.915 0.5 0.5]
% limegreen = [0.2 0.8 0.2];
% forestgreen = [0.13 0.55 0.13];
% persimmon = [0.93 0.35 0];
% OEScolor389 = [0.53, 0.0, 0.69]; %[0.53, 0.0, 0.69]; % Electric violet
% OEScolor502 = forestgreen;%[0.0 0.56 0.0];%[0.0, 0.28, 0.67];% Persian blue
% Coracolor389 = [0 0.22 0.66]; 
% Coracolor502 =  [0.83 0 0];% Azure light
% Gotocolor389 = [0.3 0.8 1]; 
% Gotocolor502 = [0.93 0.35 0];

color_init = 0.2*[1 1 1];
color_RT = [0.83 0 0];
color_OES = [0.13 0.55 0.13];

load(['OESdata_' TSfile '.mat']);
%r502 = 1e-3*r502;
%r389 = 1e-3*r389;

OESemissivity502 = f_rec502/max(f_rec502); % Normalize area to 1
% OESemissivity389 = f_rec389/trapz(r389,f_rec389');

figure; hold on; box on;
plot(r502,OESemissivity502,'Linewidth',2,'Color',color_OES)
hold on
plot(rplot,radialprofile_init,'Linewidth',1.2,'Color',color_init)
%plot(rplot,radialprofile_halo)
plot(rplot,radialprofile_RT,'Linewidth',2,'Color',color_RT)
hold off
set(gcf,'color','w')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',17) %15
xlim([-100,100])
ylim([0, 1.3])

legend('OES 502\,nm','Initial positions','RT solution','Fontsize',17,'Interpreter','latex')
xlabel('Radial position (mm)','Fontsize',21,'Interpreter','latex')
ylabel('Intensity (arb)','Fontsize',21,'Interpreter','latex')

%% Convert emission points to radial profile
% Input: radial positions of emission events such as emission_radius
% This function calculates the radial distribution of the variable
% Scaled such that peak value equals one
function profile = radial_profile(radial_positions,rplot,Nbins)
rplot = rplot*1e-3;
hist = histogram(radial_positions,Nbins); %floor(length(radial_positions)/700)
hist_radius = hist.BinEdges(1:end-1)+hist.BinWidth/2;
emissivity = hist.Values./(2*pi*hist_radius);
%emissivity = emissivity/max(emissivity); % Normalize to peak at 1
profilefct = @(r) interp1(hist_radius,emissivity,r-hist.BinWidth/2); % Return as a function -doesn't work

% Return values at rplot
profile = profilefct(rplot); 
profile(rplot>max(radial_positions)) = 0;
end