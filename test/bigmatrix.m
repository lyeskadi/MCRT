% Test big matrix nullspace solve
clear all
%close all

Nradialbins = 30;
Grr_scaling_factor = 4e6; % Temporary before actual Grr is implemented, indicates strength of reabsorption effect
Grr = Grr_scaling_factor*ones(Nradialbins);%*makeGrr(Nradialbins); %1e0*ones(Nradialbins); % TODO: replace with G from theory. Uniform G might be invalid, could cause issues with nullspace

c = crm();
c.settings.resolve_ions = 0;
c.settings.enable_ions = 0;
c.settings.enableOpacity = 0;
%c.settings.methodEIE = 'accurate';
tic;
c.run(2e12,4); % About 0.2 sec with opacity, 80msec without (Goto maybe 3 msec)
toc;
Nstates = length(c.densities);

load(['/home/lyes/Documents/PhD/Thesis/TS/matfiles/', 'TSdata_OES_TS3', '.mat']);
data=data(data.Temperature>=0.8,:); % Remove points at too low Te where CRMs do not work
nedata = data.Density*1e-6;
Tedata = data.Temperature;
radiusdata = data.Radius;

radius = linspace(radiusdata(1),radiusdata(end),Nradialbins);
nevals = interp1(radiusdata,nedata,radius);
Tevals = interp1(radiusdata,Tedata,radius);

fullmatrix = sparse(Nradialbins*Nstates,Nradialbins*Nstates);
densities_noRT = zeros(Nstates,Nradialbins);
for i = 1:Nradialbins
   c.run(nevals(i),Tevals(i));
   startind = 1+(i-1)*Nstates;
   endind = i*Nstates;
   fullmatrix(startind:endind,startind:endind)=c.RateMatrix;  
   densities_noRT(:,i) = c.densities;
end

disp('1. Check matrix without photoexcitation')
checkMatrix(fullmatrix)

tic;
%% Add photoexcitation
photoexc_matrix = sparse(Nstates,Nstates);
photoexc_matrix(5,5) = 1; % Increases 21p density proportionally to 21p density at another position
photoexc_matrix(11,11) = 1;
photoexc_matrix(19,19) = 1;

photoexc_matrix(1,5) = -1; 
photoexc_matrix(1,11) = -1;
photoexc_matrix(1,19) = -1;

disp('2. Check photoexcitation matrix')
checkMatrix(photoexc_matrix)

blocksize = Nstates;
for blockrow = 1:Nradialbins
    for blockcol = 1:Nradialbins
        block = Grr(blockrow,blockcol)*photoexc_matrix;
        
        startind_row = 1+(blockrow-1)*blocksize;
        startind_col = 1+(blockcol-1)*blocksize;

        fullmatrix(startind_row:(startind_row+blocksize-1),startind_col:(startind_col+blocksize-1)) = block + fullmatrix(startind_row:(startind_row+blocksize-1),startind_col:(startind_col+blocksize-1));
        %fullmatrix = addBlock(block, fullmatrix, blockrow, blockcol);
    end
end

disp('3. Check full matrix with photoexcitation')
checkMatrix(fullmatrix)

% Remark: if this takes long, cut time in half by using Grr' = Gr'r? or smth like this, check when G is done

%% Find the nullspace

%fulldensity = null(fullmatrix); 
%28.08: Seems ok without photoexc (Grr=0) but not 1d nullspace with it -> maybe recheck particle conservation etc 
% null() does not work with sparse matrices

% Find the convenient basis for the nullspace such that each column respresents a density vector in at a single radial position

% fulldensity_rotated = rotate_to_block_local(fulldensity, Nstates, Nradialbins);
% 
% % Check that result is indeed part of nullspace
% check_fulldensity = fullmatrix*fulldensity;
% check_fulldensity_rotated = fullmatrix*fulldensity_rotated;
% %[fulldensity,~] = eigs(fullmatrix,1,1e-12); %null(fullmatrix);
% toc;

%% Mistral

% --- Constrained Solve ---
% Construct constraint matrix: sum of densities at each position = 1 (or any fixed value)
B = zeros(Nradialbins, Nradialbins*Nstates);
for i = 1:Nradialbins
    rows = (1:Nstates) + (i-1)*Nstates;
    B(i, rows) = 1;
end
% Augmented system: [A; B] * n = [0; ones(Nradialbins,1)]
Aug = [fullmatrix; B];
b = [zeros(Nradialbins*Nstates, 1); ones(Nradialbins, 1)];
% Solve
n = Aug \ b;
% Reshape to get density vectors at each position
densities = reshape(n, Nstates, Nradialbins);

%% Plot 31P profile
plotstate = 19;

profile = densities(plotstate,:);
profile = profile/max(profile);
profile_noRT = densities_noRT(plotstate,:);
profile_noRT = profile_noRT/max(profile_noRT);

figure; hold on;
plot(radius,profile)
plot(radius,profile_noRT)
title(['Grr scaling factor = ', num2str(Grr_scaling_factor,'%.2g')])
legend([c.getname(plotstate),' With RT'], [c.getname(plotstate),' No RT'])
%% Functions

% Place block at specified position in the large matrix
% function newmat = addBlock(block,oldmat,blockrow,blockcol)
% blocksize = length(block);
% matsize = length(oldmat);
% 
% % Check if integer value
% if mod(matsize, blocksize) == 0 
%     Nblocks = matsize / blocksize;
% else
%     error('Matsize is not an integer multiple of blocksize')
% end
% 
% % Check if requested row,col is inside the large matrix
% if ~(blockrow*blocksize <= matsize && blockcol*blocksize <= matsize)
%     error('Requested row or column is outside of bounds')
% end
% 
% % Calculate block position
% startind_row = 1+(blockrow-1)*blocksize;
% startind_col = 1+(blockcol-1)*blocksize;
% 
% newmat = oldmat;
% newmat(startind_row:(startind_row+blocksize-1),startind_col:(startind_col+blocksize-1)) = block + newmat(startind_row:(startind_row+blocksize-1),startind_col:(startind_col+blocksize-1));
% end

% Make Grr matrix
function Grr = makeGrr(Nradialbins)
    Grr = zeros(Nradialbins);
    for i = 2:Nradialbins-1
        Grr(i,i)=0.4;
        Grr(i-1,i)=0.3;
        Grr(i+1,i)=0.3;
        %Grr(i-2,i)=0.1;
        %Grr(i+2,i)=0.1;        
    end
end

% Check for columns that do not sum to zero
function checkMatrix(M)
    colsums = full(sum(M,1)); % Convert to non-sparse for fprintf
    fprintf('max abs column sum = %.3e\n', max(abs(colsums)));
end

%% TODO
% 1. Resolve mfp by frequency
    % add freq roll within Doppler lineshape (write as fct)
    % calc mfp for given freq
    % When rolling distance, input mfp associated to freq
    % Reemitted freq: roll again from Doppler / Keep original one, try both and see if significant difference

% 2. Write a fct called Grr_uniform(nHe,THe) that performs the MCRT sims for a given (uniform) neutral density profile