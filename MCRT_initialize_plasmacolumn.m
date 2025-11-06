function pos = MCRT_initialize_plasmacolumn(Npoints,L_RAID,cdf_val,invcdf_val)
doPlot = 0;

%% Draw random positions distributed according to emissivity profile
% Use saved CDF^-1 from CRM_radialprofile_sampling
% if useOpacity ~= 0
%     OPtext = '_OP';
% else OPtext = []; 
% end
% load(['sampling_' filename '_' TScase OPtext '.mat']);

z = L_RAID*rand(1,Npoints); % Uniform across axial direction
r = interp1(cdf_val,invcdf_val,rand(1,Npoints)); % Use inverse sampling to reproduce radial profile expected from TS data
phi = 2*pi*rand(1,Npoints); % Cylindrical symmetry

x = r.*cos(phi);
y = r.*sin(phi);

pos = [x;y;z];

%% Optional: verify that sampling is correct - no longer compatible with current file structure ...
if doPlot
% Plot points
figure;
plot3(x,y,z,'x');
axis equal
% Compare to CRM_emissivity_radialprofile.m
[radius, profile] = CRM_emissivity_radialprofile(filename,0,TScase,useOpacity);

% Normalize the profile to PDF, account for cylindrical geometry (dV ~ 2*pi*r*dr, int(pdf)dV = 1)
profile_fun = @(r) interp1(radius,radius.*profile',r);
normalization = integral(profile_fun,0,max(radius));
normalized_profile = radius.*profile'/normalization;

figure; hold on
histogram(r,'Normalization','pdf');
plot(radius,normalized_profile); 
legend('Distribution of random samples','Expected distribution from TSdata')
end
end