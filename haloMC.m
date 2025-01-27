% Monte Carlo simulation of radiative transport
clear all
%close all

s=load('sampling.mat'); % cftval, invcftval
cftval = s.cftval;
invcftval = s.invcftval;
TScase = 'neghalf';
useOpacity = 1;

% Photon is emitted in random direction
% Travel distance is a random number drawn from the distribution f(rho)=-dT/drho (implemented by inverse sampling)
% Absorption event: roll for VUV or VIS emission. If VIS then it escapes, record how many are emitted where
% If VUV then roll another direction and emit, start again

%% Initialize variables

% Settings
TSfile = 'TS3';
filename = ['TSdata_OES_', TSfile];%'He_LCIF';%'TSdata_OES_TS2'; % Name of Thomson scattering data file for plasma bkg
Nparticles = 1e4;
R_RAID = 0.2; L_RAID = 1.5; % R_EXC=5e-2;

%nHe = 5e19; % Density of ground state atoms, assumed constant for a first attempt (replace with 20sccm value!)
% load(['f' filename]);
% ne_axis = table2array(data(1,7));
% Te_axis = table2array(data(1,9));
%c = crm(1e-6*ne_axis,Te_axis);

%ntot = ne_axis/(sum(c.poplevels(20:29))+2*c.poplevels(30));
%nground = ntot*c.poplevels(1);
nground= 1e20;
THe = 1000; % Given in Kelvin, take from LIF He data


% Fluorescence vs VUV branching ratio, neglect quenching 
Afl = 13372000;
Avuv = 566340000;
branching_ratio = Afl/(Afl+Avuv); % proba for fluorescence emission, rest is VUV
emission_pos = [];
emission_radius = [];

% Absorption strength, 1/k0 ~ mean free path
[~,k0] = Trho(11,1,nground,THe,1);

% Initial positions drawn from emissivity profile determined by TS data
initial_emission_positions = haloMC_initialize_plasmacolumn(filename,Nparticles,L_RAID,TScase,useOpacity); % Returns points according to emissivity profile predicted for TS data

% We will assume that VUV is absorbed when hitting wall
wall_absorption_pos = []; % Track where this happens

%% Main loop
parfor i = 1:Nparticles
    isVUV = 1;
    pos = initial_emission_positions(:,i)'; %initial_pos;
    dir = random_direction; % Initial direction
    while isVUV
        % Photon has been emitted, pick direction and distance 
        dx = (1/k0)*pickTravelDistance(cftval,invcftval);
        pos = pos+dir*dx; 
        if ~isInsideCylinder(pos,R_RAID,L_RAID)
            isVUV = 0; % Assume walls absorb the VUV
            wall_absorption_pos = [wall_absorption_pos;pos];
        end
        
        % Photon absorbed at pos, now roll for type of emission           
        if (isVUV & pass_test(branching_ratio))
            % Fluorescence emission
            emission_pos = [emission_pos; pos];
            emission_radius = [emission_radius sqrt(pos(1)^2+pos(2)^2)];
            %emission_radius = [emission_radius norm(pos)]; % Distance from initial position, spherical
            isVUV = 0; % Particle has finally decayed
        else
            dir = random_direction; % Another VUV is emitted from the same place
        end
    end
end

%% Tally points
N_initialemission = 1e6; %floor(branching_ratio/(1-branching_ratio)*Nparticles); % Seems to be inaccurate (was due to low number of samples)
points_initial_emission= haloMC_initialize_plasmacolumn(filename,N_initialemission,L_RAID,TScase,useOpacity); %initial_emission_positions; %generateCylinderPoints(R_EXC,L_RAID,N_initialemission);
radius_iem = sqrt(points_initial_emission(1,:).^2+points_initial_emission(2,:).^2);

rplot=1:0.01:100;

radialprofile_halo = radial_profile(emission_radius,rplot,70); % Nbins
radialprofile_halo = radialprofile_halo/max(radialprofile_halo);
radialprofile_init = radial_profile(radius_iem,rplot,100);
radialprofile_init = radialprofile_init/max(radialprofile_init);

escapeproba = length(wall_absorption_pos)/Nparticles;
%radialprofile_RT = branching_ratio*radialprofile_init + (1-branching_ratio)*radialprofile_halo; % escapeproba %@(r) branching_ratio*radialprofile_init(r) + (1-branching_ratio-escapeproba)*radialprofile_halo(r);
radialprofile_RT = branching_ratio*radialprofile_init + (1-branching_ratio-escapeproba)*radialprofile_halo;
radialprofile_RT = radialprofile_RT/max(radialprofile_RT);

if TScase == "neghalf"
    rplot =- rplot;
end

%% Plot radial profile of emission line
load(['OESdata_' TSfile '.mat']);
%r502 = 1e-3*r502;
%r389 = 1e-3*r389;

OESemissivity502 = f_rec502/max(f_rec502); % Normalize area to 1
% OESemissivity389 = f_rec389/trapz(r389,f_rec389');

figure; hold on; box on;
plot(r502,OESemissivity502)
hold on
plot(rplot,radialprofile_init)
plot(rplot,radialprofile_halo)
plot(rplot,radialprofile_RT)
hold off
set(gcf,'color','w')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',17) %15
xlim([-100,100])

legend('OES 502\,nm','No opacity','Halo','RT solution','Interpreter','latex')
xlabel('Radial position (mm)','Interpreter','latex')
ylabel('Normalized line intensity','Interpreter','latex')



%% Old
%{
% Radial distribution
% The histogram has equal size bins whereas dV = 2*pi*r*dr*L for a cylinder
% Divide bin values by 2*pi*r to get uniform dV, divide further by L and Nparticles to normalize
hist=figure;hold on;
emission_radius_all = [radius_iem, emission_radius];
h_init = histogram(radius_iem);
h_halo = histogram(emission_radius); % Halo alone
h_all = histogram(emission_radius_all); % Total = Initial + Halo
title('Radius of emission points (Divide by 2\pi r for radial density profile)')
legend('Initial emission points','Halo points','Total')
uistack(h_halo,'top')
uistack(h_all,'bottom')
% Maybe check that edge effects don't distort it, can restrict to points around TS pos

% Vertical distribution
histv=figure;hold on;
emission_verticalpos_all = [points_initial_emission(2,:), emission_pos(2,:)];
hv_init = histogram(points_initial_emission(2,:),'Normalization','pdf');
hv_halo = histogram(emission_pos(2,:),'Normalization','pdf'); % Halo alone
hv_all = histogram(emission_verticalpos_all,'Normalization','pdf'); % Total = Initial + Halo
title('Vertical position of emission points') 
legend('Initial emission points','Halo points','Total')
xlim([-0.05 0.05]);
uistack(hv_halo,'top')
uistack(hv_all,'bottom')

%% Calculate the apparent emissivity of the plasma

radius_all = h_all.BinEdges(1:end-1)+h_all.BinWidth/2;
verticalpos_all = hv_all.BinEdges(1:end-1)+hv_all.BinWidth/2;
emissivity_all = h_all.Values./(2*pi*radius_all);
emissivityv_all = hv_all.Values;

radius_halo = h_halo.BinEdges(1:end-1)+h_halo.BinWidth/2;
verticalpos_halo = hv_halo.BinEdges(1:end-1)+hv_halo.BinWidth/2;
emissivity_halo = h_halo.Values./(2*pi*radius_halo);
emissivityv_halo = hv_halo.Values;

radius_init = h_init.BinEdges(1:end-1)+h_init.BinWidth/2;
verticalpos_init = hv_init.BinEdges(1:end-1)+hv_init.BinWidth/2;
emissivity_init = h_init.Values./(2*pi*radius_init);
emissivityv_init = hv_init.Values;

%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot emission locations in 3D
figure;hold on;axis equal
%for i=1:length(emission_pos)
%plot3(0,0,0,'kx')
%plot3(wall_absorption_pos(:,1),wall_absorption_pos(:,2),wall_absorption_pos(:,3),'k+')
%plot3(pointsInside(1,:),pointsInside(2,:),pointsInside(3,:),'b+') % Initial positions
plot3(points_initial_emission(1,:),points_initial_emission(2,:),points_initial_emission(3,:),'b+') % Emission position of particles that originally decayed by fluorescence
plot3(emission_pos(:,1),emission_pos(:,2),emission_pos(:,3),'rx')
%end



%% Plot emissivity profiles
%{
figure; hold on;
plot(radius_all,emissivity_all)
plot(radius_halo,emissivity_halo)
plot(radius_init,emissivity_init)
legend('Total emissivity','Emissivity of halo','Initial emissivity')
title('Radial profile')

figure; hold on;
plot(verticalpos_all,emissivityv_all)
plot(verticalpos_halo,emissivityv_halo)
plot(verticalpos_init,emissivityv_init)
legend('Total emissivity','Emissivity of halo','Initial emissivity')
title('Vertical profile')
xlim([-0.05,0.05])

% Plot normalized emissivity profiles
profile_all = @(r) interp1(radius_all,emissivity_all,r);
profile_halo = @(r) interp1(radius_halo,emissivity_halo,r);
profile_init = @(r) interp1(radius_init,emissivity_init,r);
% normalization_all = integral(profile_all,0,max(radius));
% normalization_halo = integral(profile_halo,0,max(radius));
% normalization_init = integral(profile_init,0,max(radius));
% pdf = @(r) profile_all(r)/normalization_all;

normalized_emissivity_profile_all = normalize_profile(profile_all,[0, 0.05]);
normalized_emissivity_profile_halo = normalize_profile(profile_halo,[0, 0.05]);
normalized_emissivity_profile_init = normalize_profile(profile_init,[0, 0.05]);

rplot = linspace(0,0.1,100);
figure; hold on; box on;
plot(rplot,normalized_emissivity_profile_init(rplot))
plot(rplot,normalized_emissivity_profile_halo(rplot))
plot(rplot,normalized_emissivity_profile_all(rplot))
legend('Initial emission points','Halo')
title('Radial profile')
%}

%% Calculate CRM prediction with alternative settings
% For RT calcs start from no opacity case, but for results show: escape factor vs RT vs exp data
[radiusCRM, profileCRM] = CRM_emissivity_radialprofile(filename,0,1); % Enable opacity for calcs
%CRMemissivity = @(x) interp1(radiusCRM,profileCRM,abs(x));
%% Compare to OES data
%close all
load(['OESdata_' TSfile '.mat']);
r502 = 1e-3*r502;
r389 = 1e-3*r389;

OESemissivity502 = f_rec502/trapz(r502,f_rec502); % Normalize area to 1
OESemissivity389 = f_rec389/trapz(r389,f_rec389'); 
RTemissivity = 0.5*emissivity_all/trapz(radius_all,emissivity_all);% 2*pi*radius_all.*
INITemissivity = 0.5*emissivity_init/trapz(radius_init,emissivity_init);
CRMemissivity = 0.5*profileCRM/trapz(radiusCRM,profileCRM);

%OESemissivity502 = f_rec502/max(f_rec502); % Peak scaled to 1 as in Cora chapter
%OESemissivity389 = f_rec389/trapz(r389,f_rec389');

%radialprofile_halo = emissivity_halo/max(emissivity_halo); % Peak scaled to 1
%radialprofile_init = emissivity_init/max(emissivity_init);

% Account for lost vuv photons on the wall: they don't contribute to halo, and the initial emission is not just 2% of total emission but more
%escapeproba = length(wall_absorption_pos)/Nparticles;
%RTemissivity = branching_ratio*radialprofile_init + (1-branching_ratio-escapeproba)*radialprofile_halo;

%INITemissivity = 0.5*emissivity_init/trapz(radius_init,emissivity_init);
%CRMemissivity = 0.5*profileCRM/trapz(radiusCRM,profileCRM);

% Compare 502nm line measurement vs prediction
mirror = 1;
if mirror
    radius_init = sort([-radius_init radius_init]);
    radius_all = sort([-radius_all radius_all]);
    radiusCRM = sort([-radiusCRM' radiusCRM']);
    RTemissivity = [RTemissivity(end:-1:1), RTemissivity];
    CRMemissivity = [CRMemissivity(end:-1:1), CRMemissivity];    
end

%% Plot OES vs CRM
azure = [0 0.5 1];
%azurelight = [0.3 0.8 1];
rosso = [0.83 0 0];
%cadmium = [0.89 0 0.13];

figure; hold on; box on;
plot(r502,OESemissivity502,'LineWidth',1.2, 'Color', azure)
%plot(radius_init,CRMemissivity,'LineWidth',1.2)
plot(radiusCRM,CRMemissivity,'k-','LineWidth',1.2)
plot(radius_all,RTemissivity, 'LineWidth',1.2, 'Color', rosso)
xlabel('Radial position (m)','Interpreter','latex')
ylabel('Normalized line intensity','Interpreter','latex')
legend('OES 502nm', 'CRM with EF.','CRM with RT','Interpreter','latex')
title('502\,nm line','Interpreter','latex');
set(gcf,'color','w')
%set(h,'LineWidth',1.2)
%set(h2,'LineWidth',1.2)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',17) %15

% Check 389nm line
% [crm_radius389, crm_profile389] = CRM_emissivity_radialprofile(filename,0);
% 
% crm_r389 = linspace(crm_radius389(1),crm_radius389(end),37);
% crm_pr389 = interp1(crm_radius389,crm_profile389,crm_r389);
% crm_pr389 = 0.5*crm_pr389/trapz(crm_r389,crm_pr389);
% figure; hold on; box on;
% plot(r389,OESemissivity389)
% plot(crm_r389,crm_pr389);
% xlabel('Radial position (m)')
% ylabel('Line intensity')
% legend('OES 389nm', 'CRM prediction')

%% Fit exponential to halo
%{
expfit_bounds = [0.04 max(r502)]; 
expfit_inds = find(radius_all>expfit_bounds(1) & radius_all<expfit_bounds(2));
expfit_radius = radius_all(expfit_inds);
expfit_RTemissivity = RTemissivity(expfit_inds);
expfit_OESemissivity = interp1(r502,OESemissivity502,expfit_radius);

% params(1) = amplitude, params(2) = exp decay length
fitfct = @(params,x) params(1)*exp(-(x-expfit_bounds(1))/params(2));
startParams = [max(expfit_OESemissivity/2), 1/k0];
paramsRT = lsqcurvefit(fitfct, startParams, expfit_radius, expfit_RTemissivity);
paramsOES = lsqcurvefit(fitfct, startParams, expfit_radius, expfit_OESemissivity);

fitplotr = linspace(expfit_radius(1),expfit_radius(end),100);

figure; hold on; box on;
plot(r502,OESemissivity502)
%plot(r389,OESemissivity389)
%plot(radius_init,CRMemissivity)
plot(radius_all,RTemissivity)
plot(fitplotr,fitfct(paramsRT,fitplotr))
plot(fitplotr,fitfct(paramsOES,fitplotr))
xlabel('Radial position (m)')
ylabel('Line intensity')
legend('OES 502nm','CRM predicition with RT',['Exponential fit to RT, L=' num2str(paramsRT(2)*k0) ' 1/k0'],['Exponential fit to OES, L=' num2str(paramsOES(2)*k0) ' 1/k0'])
%}

%}
%% Function definitions
function rdir = random_direction()
    randvect = rand(1,3)-0.5*[1,1,1];
    rdir = randvect/norm(randvect);
end

function didPass = pass_test(proba_pass)
    roll = rand;
    if roll < proba_pass
        didPass = 1;
    else
        didPass = 0;
    end
end

function  res = dTdrho(k0,rho)
    for i = 1:length(rho)
        integrand = @(x) exp(-2*x.^2-k0*rho(i)*exp(-x.^2));
        res(i) = -(k0/sqrt(pi))*integral(integrand,-inf,inf);
    end
end

function k0rho = pickTravelDistance(cftval,invcftval)
    uniform_sample = rand();
    k0rho = interp1(cftval,invcftval,uniform_sample);
end

function res = isInsideCylinder(pos,r,L)
if(sqrt(pos(1)^2+pos(2)^2) > r || pos(3) > L || pos(3) < 0)
    res = 0;
else
    res = 1;
end
end

function normd = normalize_profile(profile,lims)
normalization = integral(profile,lims(1),lims(2));
normd = @(r) profile(r)/normalization;
end

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