%%% A MATLAB CODE FOR THE FILTER- BASED LEVEL SET TOPOLOGY METHOD
function filter_based_levelset(nelx, nely,V_max, penal,rmin) %% PARAMETERSS
E0 = 1;      % Young modulus assigned to solid domain 
ep_E = 1e-4; % Mimicing parameter for penalizing Young modulus over the void domain 
nu=0.3;      % Poisson ratio 
PI=[10, 8];  % [Proportional gain, Integral gain] 
I_E=0;       % Initial value for the integral part of Eq 
dt=0.1;      % Fictitious time step 
T=6;         % The time at which volume fraction is expected to reach V max 
Tinf=20;     % The total time for optimization process. 
             % The total iterations is obtained from Tinf/dt (see line 36 of the code)
V_inp=max(1-(1-V_max)/T* [dt:dt:Tinf], V_max); % Linear input reference (Yaghmaei2019) 
phi_n_v=ones((nely+1)*(nelx+1),1)*.5;  %Initial level set function for each node in vector form 
% FILTER PREPARATION SEE Andreassen et al. (2011)
[dy, dx] = meshgrid(-ceil(rmin)+1: ceil (rmin)-1, -ceil(rmin)+1: ceil(rmin)-1);
H= max(0,rmin-sqrt(dx.^2+dy.^2))/sum(sum(max(0,rmin-sqrt(dx.^2+dy.^2))));  %filter kernel 
% FINITE ELEMENT ANALYSIS SEE Andreassen et al. (2011)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE0 = E0/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
KME0=E0/(1-nu^2)^2/24*((3*nu^2-2*nu+3)*[A11 A12;A12' A11]+(-nu^2+6*nu-1)*[B11 B12;B12' B11]); 
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
edofVec = reshape(2*nodenrs (1: end-1, 1: end-1)+1, nelx*nely, 1);
edofMat = repmat (edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1);
iK = reshape(kron(edofMat, ones(8,1))',64*nelx*nely, 1);
jK = reshape(kron(edofMat, ones(1,8))',64*nelx*nely, 1);
% DEFINE LOADS AND SUPPORTS (CANTILEVER BEAM)
F= sparse([2*((nely+1)*(nelx+1)-nely/2)],[1],[1],2*(nely+1)*(nelx+1),1);
U = zeros(2* (nely+1)* (nelx+1), 1);
fixeddofs=1: 2*(nely+1);
alldofs = [1:2* (nely+1) * (nelx+1) ]; 
freedofs = setdiff(alldofs, fixeddofs);
%% OPTIMIZATION LOOP
for It=1: Tinf/dt
    % FILTERING OF THE LEVEL SET FUNCTION AND UPDATING MATERIAL DISTRIBUTION
    phi_n_m=reshape(phi_n_v, nely+1, nelx+1); % Level set function for each node in matrix form 
    phi_n_m=imfilter(phi_n_m, H,'same');       % Filtering of the level set function. Refer to "Matlab Online Help, [online] Available: http: //www.mathworks. com Initial" for more details on the function imfilter.
    phi_e_m=0.25*(phi_n_m(1: end-1,1: end-1)+phi_n_m (2: end,1: end-1) +phi_n_m(1: end-1,2: end) +phi_n_m (2: end, 2: end) ) ; % Mapping the level set function onto elements in matrix form (Otomori 2015).
    str(:,:) = sigmf(phi_e_m(:,:), [penal+It*0 0]); % Updating material distribution. Refer to "Matlab Online Help,[online]Available:http://www. mathworks.comInitial"formore details on the functionsigmf.
    V = sum(str(: ))/(nely*nelx);
    % FE-ANALYSIS, OBJECTIVE FUNCTION AND TOPOLOGICAL DERIVATIVE 
    sK= reshape(KE0(:)* (ep_E +str(:)'*(1-ep_E)),64*nelx*nely,1) ;
    K=sparse (iK, jK, sK); K=(K+K')/2; % See Andreassen et al. (2011) 
    U(freedofs) = K(freedofs, freedofs)\F(freedofs);
    J1(It)=full(F(2*((nely+1)*(nelx+1)-nely/2),1))*U(2*((nely+1)*(nelx+1)-nely/2),1);  % Calculating the objective function
    alfa = 1/J1 (It) * (nely+1) * (nelx+1) ;
    TDJ_e=(ep_E+str*(1- ep_E)).*reshape(sum((U(edofMat)*KME0).*U(edofMat), 2), nely, nelx); % Calculating the topological derivative on each element.
    td2=[TDJ_e(1,1) TDJ_e(1,:) TDJ_e(1,end); TDJ_e(:,1) TDJ_e TDJ_e(:,end); TDJ_e(end,1) TDJ_e(end,:) TDJ_e(end, end)]; 
    % UPDATING THE LAGRANGE MULTIPLIER
    E = (V - V_inp(It)) ; 
    I_E = I_E+ E*dt;
    Lambda = PI(1) * E * min(1+It/10,3)+ PI(2) * I_E; 
    % UPDATING THE LEVEL SET FUNCTION AND PRINTING RESULTS
    TDJ_n=0.25*(td2(1:end-1,1:end-1) +td2(2:end, 1:end-1) +td2(1:end-1,2:end) +td2(2:end,2:end)); % Mapping the topological derivative onto nodes in matrix form (Otomori 2015).
    phi_n_v=phi_n_v+(TDJ_n(:) *alfa-Lambda) *dt; %Updating the level set function in the vector form.
    phi_n_v=min(1, max(-1, phi_n_v));
    disp(['t.:' sprintf('%2.1f',It*dt) ' Jc/Jc_0:' sprintf('%2.3f',J1(It)/2)...
          ' Vol.:' sprintf('%2.2f',V/1)])
    colormap(gray); imagesc(1-(phi_e_m>0 | str==1)); caxis([0 1]); axis equal; axis off; drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @ARTICLE{Yaghmaei2020,
%   author = {Mohammad Yaghmaei and Ali Ghoddosian and Mohammad Mahdi Khatibi},
%   title = {A filter-based level set topology optimization method using 
%            a 62-line MATLAB code},
%   journal = {Structural and Multidisciplinary Optimization},
%   year = {2020},
%   volume = {62},
%   pages = {1001-1018}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%