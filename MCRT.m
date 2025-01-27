% Monte Carlo radiative transfer simulation
function [emission_radius, escapeproba] = MCRT(nground,branching_ratio,TSfile, Nparticles,TScase,useOpacity)
% Load sampling of T(rho) CDF^-1, used to generate the distance a photon travels before absorption
s=load('sampling.mat'); % cftval, invcftval
cftval = s.cftval;
invcftval = s.invcftval;
%TScase = 'half';
%useOpacity = 1;

% Photon is emitted in random direction
% Travel distance is a random number drawn from the distribution f(rho)=-dT/drho (implemented by inverse sampling)
% Absorption event: roll for VUV or VIS emission. If VIS then it escapes, record how many are emitted where
% If VUV then roll another direction and emit, start again

%% Initialize variables

% Settings
%TSfile = 'TS3';
filename = ['TSdata_OES_', TSfile];%'He_LCIF';%'TSdata_OES_TS2'; % Name of Thomson scattering data file for plasma bkg
%Nparticles = 1e6;
R_RAID = 0.2; L_RAID = 1.5; % R_EXC=5e-2;

%nHe = 5e19; % Density of ground state atoms, assumed constant for a first attempt (replace with 20sccm value!)
% load(['f' filename]);
% ne_axis = table2array(data(1,7));
% Te_axis = table2array(data(1,9));
%c = crm(1e-6*ne_axis,Te_axis);

%ntot = ne_axis/(sum(c.poplevels(20:29))+2*c.poplevels(30));
%nground = ntot*c.poplevels(1);
%nground= 1e20;


THe = 1000; % Given in Kelvin, take from LIF He data

% Fluorescence vs VUV branching ratio, neglect quenching 
 % proba for fluorescence emission, rest is VUV
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
    
    % Include initial emission (rewrite later)
    if pass_test(branching_ratio)
        emission_pos = [emission_pos; pos];
        emission_radius = [emission_radius sqrt(pos(1)^2+pos(2)^2)];
        isVUV = 0;
    end
    
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

%% Set outputs
escapeproba = length(wall_absorption_pos)/Nparticles;
end

%% Function definitions

% Roll a random direction for the test particle
function rdir = random_direction()
    randvect = rand(1,3)-0.5*[1,1,1];
    rdir = randvect/norm(randvect);
end

% Roll a binary choice with proba p and 1-p
function didPass = pass_test(proba_pass)
    roll = rand;
    if roll < proba_pass
        didPass = 1;
    else
        didPass = 0;
    end
end

% Roll the distance traveled by the photon before absorption
function k0rho = pickTravelDistance(cftval,invcftval)
    uniform_sample = rand();
    k0rho = interp1(cftval,invcftval,uniform_sample);
end

% Check if particle is inside RAID
function res = isInsideCylinder(pos,r,L)
if(sqrt(pos(1)^2+pos(2)^2) > r || pos(3) > L || pos(3) < 0)
    res = 0;
else
    res = 1;
end
end

