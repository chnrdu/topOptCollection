% Matlab code for topology optimization using a reaction diffusion equation
function [str,phi] = levelsetRD(nelx,nely,Vmax,tau)
%% Parameter definition
E0 = 1;
Emin = 1e-4;
nu = 0.3;
nvol = 100;
dt = 0.1;
d = -0.02;
p = 4;
phi = ones((nely+1)*(nelx+1),1);
str = ones(nely,nelx);
volfrac = sum(str)/(nelx*nely);
%% Finite element analysis preparation
% For displacement field
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 3; 9 4 3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
C = 1/((1-nu^2)/24)*([A11 A12;A12' A11] + nu*[B11 B12;B12' B11]);
a1 = 3*(1-nu)/(2*(1+nu)*(7-5*nu))+(-(-14*nu+15*nu^2)*E0)/(1-2*nu)^2;
a2 = -6*(1-nu)/(2*(1+nu)*(7-5*nu))*E0;
A = (a1*a2)/24*([A11 A12;A12' A11]+(a1/(a1*a2))*[B11 B12;B12' B11]);
nodeMat = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodeMat(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2 -1 -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
MMdif = 1/6*[ 4 -1 -2 -1;-1 4 -1 -2;-2 -1 4 -1;-1 -2 -1 4];
Dif = 1/36*[ 4 2 1 2;2 4 2 1; 1 2 4 2; 2 1 2 4];
edofVec2= reshape(nodeMat(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat2= repmat(edofVec2,1,4)+repmat([0 nely+1 nely -1],nelx*nely,1);
iM = reshape(kron(edofMat2,ones(4,1))',16*nelx*nely,1);
jM = reshape(kron(edofMat2,ones(1,4))',16*nelx*nely,1);
sM = reshape(kron(ones(1,nelx*nely),ones(1,4*4)),16*nelx*nely,1);
MM = sparse(iM,jM,sM);
sMMdif = reshape(kron(ones(1,nelx*nely),reshape(MMdif,1,16)),16*nelx*nely,1);
MMdif = sparse(iM,jM,sMMdif);
%% Loads and boundary settings
F = sparse(2*(nely+1)*(nelx+1),1);
K = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
F((nely+1)*(nelx+1)*2*(-round(nely/32)+1):2:(nely+1)*(nelx+1)*2*(-round(nely/32)+1),1) = 1;
fixeddofs = 1:2*(nely+1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
fixeddofs_phi = [];
freedofs_phi = sort([1:nely+1 nely+2:nely+1:(nely+1)*(nelx) 2*(nely+1):nely+1:(nely+1)*(nelx) ...
    (nely+1)*nelx+1:(nely+1)*(nelx+1)]);
phi(fixeddofs_phi) = 0;
alldofs_phi = 1:(nely+1)*(nelx+1);
freedofs_phi = setdiff(alldofs_phi,fixeddofs_phi);
%% Main loop
for iterNum = 1:1000
    %% FE-analysis, calculate sensitivities
    sK = reshape((Emin+(E0-Emin)*str(:)'.^p),64*nelx*nely,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs) \ F(freedofs);
    U(fixeddofs) = 0;
    se = (Emin+(str(:)'.^p)*(E0-Emin)).*reshape(sum(reshape(U(edofMat),nely,nelx) ...
    .^2),nely*nelx,1);
    DU = (1e-8+str(:)'.^(1-1e-8)).*reshape(sum(U(edofMat).^2),nely*nelx,1);
    OBJ = (-0.5*(DU(:)').^1)*se;
    phi(freedofs_phi) = (MM(freedofs_phi,freedofs_phi)+tau*dt*MMdif(freedofs_phi,freedofs_phi)) ...
        \ (MM(freedofs_phi,freedofs_phi)*phi(freedofs_phi)- ...
        dt*d*MM(freedofs_phi,freedofs_phi)*ones(length(freedofs_phi),1));
    objective(iterNum) = sum(OBJ(:));
    vol = sum(str(:))/(nelx*nely);
    %% Print results
    disp(['It.: ' num2str(iterNum) ' Compl.: ' sprintf('%10.4e',objective(iterNum)/(nelx*nely))...
        ' Vol.: ' sprintf('%6.2f',vol)]);
    colormap(gray); imagesc(-str,[-1,0]); axis equal; axis tight; axis off; drawnow;
    %% Convergence check
    if iterNum>nvol && (abs(vol-Vmax)<0.005) && all(abs(objective(end-3:end-1) - objective(end)) < 0.01*abs(objective(end)))
        return;
    end
    if iterNum>nvol && (abs(vol-Vmax)<0.005) || (iterNum==1000)
        return;
    end
end

%%----------------------------------------------------------------------------
% Title: Matlab code for a level set-based topology optimization method 
%        using a reaction diffusion equation
% authors: Masaki Otomori, Takayuki Yamada, Kazuhiro Izui & Shinji Nishiwaki 
% journal: Structural and Multidisciplinary Optimization
% Volume 51, pages 1159–1172, (2015)
% https://doi.org/10.1007/s00158-014-1190-z
%%----------------------------------------------------------------------------

