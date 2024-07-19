% Proportional Topology Optimization compliance (PTOc) - Half MBB beam - (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = PTOc_mbb(alpha,E0,Emin,L,lv,ld,nelx,nely,nu,penal,rmin,vlim,xlim)
% Setup Finite Element Analysis
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% Define Loads and Supports
iF = 2*(nely+1)*(0:ld-1)+2; 
jF = ones(1,ld); 
sF = -lv/ld*ones(ld,1);
F = sparse(iF,jF,sF,2*(nely+1)*(nelx+1),1);
% Define Displacement and DOF Sets
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1),2*((nelx+1)*(nely+1)-ld+1:(nelx+1)*(nely+1)));
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
% Setup Stress Analysis
B = (1/2/L)*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1; -1 -1 -1 1 1 1 1 -1];
DE = (1/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
% Setup Filter
iW = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jW = ones(size(iW));
sW = zeros(size(iW));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iW(k) = e1;
        jW(k) = e2;
        sW(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
w = sparse(iW,jW,sW);
W = bsxfun(@rdivide,w,sum(w,2));
% Initialize Iteration
x = repmat(vlim,nely,nelx);
xNew = zeros(size(x));
loop = 0;
change = Inf;
% Run Iteration
while (true)
 loop = loop+1;
 % Finite Element Analysis
 E = Emin+x(:)'.^penal*(E0-Emin);
 sK = reshape(KE(:)*E,64*nelx*nely,1); 
 K = sparse(iK,jK,sK); K = (K+K')/2;
 U(freedofs) = K(freedofs,freedofs)\F(freedofs);
 % Stress Calculation
 s = (U(edofMat)*(DE*B)').*repmat(E',1,3);
 vms = reshape(sqrt(sum(s.^2,2)-s(:,1).*s(:,2)+2.*s(:,3).^2),nely,nelx);
 % Compliance Calculation 
 ce = E'.*sum((U(edofMat)*KE).*U(edofMat),2);    
 C = reshape(ce,nely,nelx);
 % Print Results
 fprintf('It:%5i Max_vms:%5.2f Comp:%8.2f Vol:%5.2f Ch:%6.3f\n',...
         loop,max(vms(:)),sum(C(:)),mean(x(:)),change);
 % Plot Results
 colormap(flipud(gray)); 
 subplot(2,1,1); imagesc(x); axis equal off; text(2,-2,'x');
 subplot(2,1,2); imagesc(C); axis equal off; text(2,-2,'C'); drawnow;
 % Check Stop Criteria
 if(change < 0.01 && loop > 50); break; end;
 % Optimization Algorithm (PTOc)
 xTarget = nelx*nely*vlim;
 xRemaining = xTarget;
 xNew(:) = 0;
 C_proportion = C/sum(C(:));
 while (xRemaining > 0.001)
  xDist = xRemaining.*C_proportion;
  xNew(:) = xNew(:)+W*xDist(:);
  xNew = max(min(xNew,xlim(2)),xlim(1));
  xRemaining = xTarget-sum(xNew(:));
 end
 x = alpha*x+(1-alpha)*xNew;
 change = max(abs((1/alpha-1)*(xNew(:)-x(:))));
end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015 University of Pittsburgh. All rights reserved.
% 
% Any person who obtained a copy of this software can (in part or whole) copy, 
% modify, merge, publish, and distribute the software on condition of retaining
% this license with the software. The user is allowed to utilize the software
% for all purposes but commercial. Also, appropriate credit must be provided.
% 
% The software is provided "as is", without warranty of any kind, express or 
% implied, including but not limited to the warranties of merchantability, 
% fitness for a particular purpose and noninfringement. In no event shall the 
% authors or copyright holders be liable for any claim, damage or other 
% liability, whether in an action of contract, tort or otherwise, arising from, 
% out of or in connection with the software or the use or other dealing in the 
% software.
% 
% The software is coded by Emre Biyikli (biyikli.emre@gmail.com) and Albert C. 
% To (albertto@pitt.edu). The software is substantially inherited from 
% Andreassen E, et al. "Efficient topology optimization in MATLAB using 88 lines
% of code." Structural and Multidisciplinary Optimization 43.1 (2011): 1-16. 
% The software can be downloaded from www.ptomethod.org. 
% The journal article to the software is ...
% 

