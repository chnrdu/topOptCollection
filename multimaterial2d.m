%% LOOP PARAMETERS
clear;
NumIter = 300;delta=1e-9;
nelx=200; nely=100;
penal=3; rmin=5; NMaterial=3;
NCON=NMaterial; NDV=nelx*nely*NMaterial; NDV3=nelx*nely;
movelimit=0.1; beta = 1;
%% MATERIAL PROPERTIES
par.E(1)=1; par.E(2)=2; par.E(3)=5;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
lk = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
loopbeta=0;
x=ones(NDV,1)*0.5;
for i=0:NMaterial-1
  xTilde(1:nely,1:nelx,i+1)=reshape(x(1+i*NDV3:NDV3+i*NDV3),nely,nelx);
end
xPhys = (tanh(beta*0.5) + tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5) + tanh(beta*(1-0.5)));
%% MMA PARAMETER INITIALIZE
a0=0; a=zeros(NCON,1);cc=1.0e6*ones(NCON,1); d=0*ones(NCON,1);
xmin=1e-3*ones(NDV,1); xmax=1*ones(NDV,1);
xold=x; xolder=xold;low=0; upp=1;
%% START ITERATION
for iter=1:NumIter
  loopbeta=loopbeta+1;
  %% FE-ANALYSIS
  [KE0] = lk; p=6;
  for i=1:NMaterial
    KE(:,:,i)=par.E(i)*KE0;
  end
  KE_test=zeros(64,nely*nelx); g=zeros(NMaterial,1);
  r1=xPhys(:,:,1); r2=xPhys(:,:,2); r3=xPhys(:,:,3);
  for j=1:NMaterial
    post_x(1:nely,1:nelx,j)=(xPhys(:,:,j).*(r1.^p + r2.^p + r3.^p ).^(1/p))./(r1 + r2 + r3 + delta);
    dr(1:nely,1:nelx,j)=(post_x(1:nely,1:nelx,j)).^penal;
    KE_test=reshape(KE(:,:,j),64,1)*reshape(dr(1:nely,1:nelx,j),1,nely*nelx)+KE_test;
    g(j)=sum((sum(xPhys(:,:,j))'));
  end
  sK = reshape(KE_test,64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  %% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
  F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
  F(2,1) = -0.1;
  fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
  alldofs = [1:2*(nely+1)*(nelx+1)];
  freedofs = setdiff(alldofs,fixeddofs);
  %% SOLVING
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  U(fixeddofs,:)= 0;
  c=F'*U;
  %% SENSITIVITY ANALYSIS
  dc=zeros(nely,nelx,NMaterial);
  dgdx_test=zeros(NMaterial,nely*nelx*NMaterial);
  KES=reshape(sum(U(edofMat)*KE0.*U(edofMat),2),nely,nelx);
  for n=1:NMaterial
    for m=1:NMaterial
      if m==n
        dmdrn(1:nely,1:nelx,m)=penal.*post_x(1:nely,1:nelx,m).^(penal -1).*...
          ((r1.^p + r2.^p + r3.^p).^(1/p)./(r1 + r2 + r3 + delta) -(xPhys(:,:,m).*(r1.^p + r2.^p + r3.^p).^(1/p))./(r1 + r2 + r3 + delta).^2 +...
          (xPhys(:,:,m).*xPhys(:,:,n).^(p -1).*(r1.^p + r2.^p + r3.^p).^(1/p -1))./(r1 + r2 + r3 + delta));
        mass_dmdrn=1;
      else
        dmdrn(1:nely,1:nelx,m)=penal.*post_x(1:nely,1:nelx,m).^(penal -1).*...
          (((xPhys(:,:,m).*xPhys(:,:,n).^(p -1).*(r1.^p + r2.^p + r3.^p).^(1/p -1))./(r1 + r2 + r3 + delta)- (xPhys(:,:,m).*(r1.^p + r2.^p + r3.^p).^(1/p))./(r1 + r2 + r3 + delta).^2));
        mass_dmdrn=0;
      end
      dgdx_test(m,1+(n-1)*NDV3:NDV3+NDV3*(n-1))=repmat(mass_dmdrn,1,nelx*nely);
    end
    dc(:,:,n)=-dmdrn(:,:,1).*KES.*par.E(1)-dmdrn(:,:,2).*KES.*par.E(2)-dmdrn(:,:,3).*KES.*par.E(3);
  end
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  filtereddc=[];dfdx=[];
  for i=1:NMaterial
    dxx = beta* (1-tanh(beta*(xTilde-0.5)).*tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5) + tanh(beta*(1-0.5)));
    [filtereddc(:,:,i)] = H*(reshape(dc(:,:,i).*dxx(:,:,i),NDV3,1)./Hs);
    dfdx=[dfdx;filtereddc(:,:,i)];
    g(i)=(g(i)/NDV3-0.5/3)*1;
    dgdx1=(dgdx_test(i,:)/NDV3)*1;
    for j=0:NMaterial-1
      dgdx2=reshape(dgdx1(:,1+j*NDV3:NDV3+j*NDV3),nely,nelx);
      dgdx(i,1+j*NDV3:NDV3+j*NDV3) = H*(reshape(dgdx2(:,:).*dxx(:,:,i),NDV3,1)./Hs);
    end
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  f=c*1e2; dfdx=dfdx*1e2; cons(iter,1)=g(1);cons(iter,2)=g(2);cons(iter,3)=g(3);
  for j=1:length(x)
    xmin(j)=max(1e-3,x(j)-movelimit); xmax(j)=min(1,x(j)+movelimit);
  end
   %% MMA OPTIZATION
  [xnew, y, z, lamda, ksi, eta, mu, zeta, s, low, upp] = ...
    mmasub(NCON,length(x), iter, ...
    x, xmin, xmax, xold, xolder, ...
    f, (dfdx), g, dgdx, low, upp, a0, a, cc, d);
  %% PRINT RESULTS
  disp(sprintf('Iter.:%3d Obj.: %8.4e max constraint: %6.3e', iter, f,max(g)))
  %% UPDATE PARAMETER
  change=norm(abs(xnew-x));
  xolder=xold; xold=x; x = xnew;
  for i=0:NMaterial-1
    xTilde(1:nely,1:nelx,i+1) = reshape((H*xnew(1+i*NDV3:NDV3+i*NDV3,:))./Hs,nely,nelx);
  end
  xPhys = (tanh(beta*0.5) + tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5) + tanh(beta*(1-0.5)));
  if beta < 100 && (loopbeta >= 50)
    beta = 2*beta;loopbeta = 0;
    fprintf('Parameter beta increased to %g.\n',beta);
  end
end
%% PLOT DENSITIES
xx=post_x(:,:,1);
for i=2:NMaterial
  xx=xx+i*post_x(:,:,i);
end
imagesc(xx); axis equal; axis off;colorbar; hold on;
plot([0 nelx nelx 0 0]+0.5,[0 0 nely nely 0]+0.5,'k');axis([-1 nelx+1 -1 nely+1])
saveas(figure(1),'result.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zheng R et al. (2024), Applied Sciences, DOI: 10.3390/app14020657           %
% Title:                                                                      %
%   An Efficient Code for the Multi-Material Topology Optimization            %
%   of 2D/3D Continuum Structures Written in Matlab                           %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
