% Calculate Grr matrix
clear all
close all
%% 0. Settings and load required data
troubleshooting = 0; % Enable if you want to check out some plots
save_data = 1;

nHe = 7e19;
THe = 1000;

absorbing_state = 1;
[~,k0_21p] = Trho(5,absorbing_state,nHe,THe,1);
[~,k0_31p] = Trho(11,absorbing_state,nHe,THe,1);
[~,k0_41p] = Trho(19,absorbing_state,nHe,THe,1);

s=load('sampling.mat'); % cftval, invcftval
cftval = s.cftval;
invcftval = s.invcftval;

%% 1. Generate points in RAID cylinder
Nradialbins = 100;
Npoints = 1e7; % Number of points in each bin

L_RAID = 1.5;
R_RAID = 0.2;

% z = L_RAID*rand(1,Npoints); % Uniform across axial direction
% r = R_RAID*rand(1,Npoints); % Same number of points in each shell -> lower density at larger radii
% phi = 2*pi*rand(1,Npoints); % Cylindrical symmetry
% 
% x = r.*cos(phi);
% y = r.*sin(phi);
% 
% pos = [x;y;z];


dr = R_RAID/Nradialbins;
r_edges = 0:dr:R_RAID;
bin_volumes = pi*L_RAID*((r_edges(1:end-1)+dr).^2-r_edges(1:end-1).^2);%2*pi*L_RAID*dr*(dr/2+r_edges(1:end-1));
Grr_21p = zeros(Nradialbins);
Grr_31p = zeros(Nradialbins);
Grr_41p = zeros(Nradialbins);

rdir = random_directions(Npoints); % Precompute random directions once, since this is surprisingly costly to do
% Further speed improvements: move pickTravelDistance calls here as well and only calculate once instead of Nradialbins times!
% However, this requires being careful with possible correlations. Example: break/reduce correlations by permutating the directions/distances inside the loop
% => Worth doing this at some point, but requires testing

parfor i = 1:Nradialbins
%for i = 1:Nradialbins
    % Generate initial positions
    z_init = L_RAID*rand(1,Npoints);
    r_init = (i-1)*dr + dr*rand(1,Npoints);
    phi_init = 2*pi*rand(1,Npoints);
    x_init = r_init.*cos(phi_init);
    y_init = r_init.*sin(phi_init);
    pos_init = [x_init;y_init;z_init]';
    
    % Generate displacement vectors
    rho_21p = pickTravelDistance(Npoints,cftval,invcftval,k0_21p);
    rho_31p = pickTravelDistance(Npoints,cftval,invcftval,k0_31p);
    rho_41p = pickTravelDistance(Npoints,cftval,invcftval,k0_41p);
    
    % rdir = random_directions(Npoints); -> Moved outside loop for speed, although could cause issues due to correlations!
    displacement_21p = rdir.*rho_21p';
    displacement_31p = rdir.*rho_31p';
    displacement_41p = rdir.*rho_41p';
    
    % Calculate final positions
    pos_final_21p = pos_init + displacement_21p;
    pos_final_31p = pos_init + displacement_31p; 
    pos_final_41p = pos_init + displacement_41p;
    
    % Remove points that left axially
    pos_final_21p((pos_final_21p(:,3)>L_RAID | pos_final_21p(:,3)<0),:) = [];
    pos_final_31p((pos_final_31p(:,3)>L_RAID | pos_final_31p(:,3)<0),:) = [];
    pos_final_41p((pos_final_41p(:,3)>L_RAID | pos_final_41p(:,3)<0),:) = [];
    
    % Calcuate radius
    r_final_21p = sqrt(pos_final_21p(:,1).^2 + pos_final_21p(:,2).^2);
    r_final_31p = sqrt(pos_final_31p(:,1).^2 + pos_final_31p(:,2).^2);
    r_final_41p = sqrt(pos_final_41p(:,1).^2 + pos_final_41p(:,2).^2);
    
    % Find which bin each point belongs to. Histcounts automatically discards points at r > R_RAID
    bincounts_21p = histcounts(r_final_21p,r_edges);
    bincounts_31p = histcounts(r_final_31p,r_edges);
    bincounts_41p = histcounts(r_final_41p,r_edges);
    
    % Calculate proba density, while accounting for volume difference of bins. For sanity check: look at bincounts/Npoints, must be <= 1. Expectation: close to 1 in the center, <<1 at edge
    Grr_21p(:,i) = (bincounts_21p/Npoints).*(bin_volumes(i)./bin_volumes);
    Grr_31p(:,i) = (bincounts_31p/Npoints).*(bin_volumes(i)./bin_volumes);
    Grr_41p(:,i) = (bincounts_41p/Npoints).*(bin_volumes(i)./bin_volumes);
    
    % Troubleshooting
    if troubleshooting % Change parfor into for otherwise no plot
        if i == 1
            figure;
            plot3(x_init,y_init,z_init,'bx')
            axis equal
            hold on
            plot3(pos_final_21p(:,1),pos_final_21p(:,2),pos_final_21p(:,3),'rx')
            legend('Initial position','Final position')
            xlabel('x')
            ylabel('y')
            zlabel('z')

            figure;
            plot(bincounts_21p)
            xlabel('Bin')
            ylabel('Counts')
        end
    end
end

%% Save Grr matrix and settings in a .mat file

if(save_data)
    save_name = ['Grr_data_' 'nHe' num2str(nHe,'%.e') '_Nbins_' num2str(Nradialbins) '_Npoints' num2str(Npoints,'%.e')];% 'noVolumescaling'
    save(save_name,'Nradialbins','Npoints','nHe','THe','Grr_21p','Grr_31p','Grr_41p');
end

%% Function definitions

% Roll a random direction for the test particle
% function rdir = random_direction()
%     randvect = rand(1,3)-0.5*[1,1,1];
%     rdir = randvect/norm(randvect);
% end
function rdir = random_directions(N)
    % Generate N random unit vectors in 3D
    randvect = rand(N, 3) - 0.5; % N x 3 matrix of random vectors
    norms = sqrt(randvect(:,1).^2 + randvect(:,2).^2 + randvect(:,3).^2);
    rdir = randvect ./ norms; % Normalize each row
end

% Roll the distance traveled by the photon before absorption
function rho = pickTravelDistance(N,cftval,invcftval,k0)
    uniform_sample = rand(1,N);
    k0rho = interp1(cftval,invcftval,uniform_sample);
    rho = k0rho/k0;
end

% Check if particle is inside RAID
% function res = isInsideCylinder(pos,r,L)
% if(sqrt(pos(1)^2+pos(2)^2) > r || pos(3) > L || pos(3) < 0)
%     res = 0;
% else
%     res = 1;
% end
% end