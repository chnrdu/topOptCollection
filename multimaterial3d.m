%% LOOP PARAMETERS
clear;
NumIter = 300;delta=1e-9;
nelx=60; nely=30; nelz = 15;
penal=3; rmin=5; NMaterial=3;
NDV=nelx*nely*nelz*NMaterial; NDV3=nelx*nely*nelz;
movelimit=0.1; beta =1;
%% MATERIAL PROPERTIES
par.E(1)=1; par.E(2)=2; par.E(3)=5;
nu = 0.3;
%% USER-DEFINED LOAD DOFs
il = nelx; jl = 0; kl = 0:nelz;
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);
loaddof = 3*loadnid(:) - 1;
%% USER-DEFINED SUPPORT FIXED DOFs
[jf,kf] = meshgrid(1:nely+1,1:nelz+1);
fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;
fixeddofs = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2];
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
nele = nelx*nely*nelz;
ndof=3*(nelx+1)*(nely+1)*(nelz+1);
F=sparse(loaddof,1,-1,ndof,1);
U=zeros(ndof,1);
freedofs=setdiff(1:ndof,fixeddofs);
%% PREPARE FINITE ELEMENT ANALYSIS
ke=1/(1+nu)/(2*nu-1)/144 *([-32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;3;6;8;3;3;4;-3;-3;-32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;...
	3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6; -4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;...
	-3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;...
	3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;...
	-32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;-6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;...
	-32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4; -3;-32;-6;-3;-4;-32;6;6;-32;-6;-32]+nu*[48;0;0;0;-24;-24;-12;0;-12;0;...
	24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;0;0;24;12;-12;12;0;-12;0;-12;-12;0;48;24;0;0;12;12;-12;0;24;0;-24;-24;0;...
	0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24; -12;-12;-12;-12;0;0;48;0;24;0;-24;0;-12;-12;-12;-12;12;0;0;24;12;-12;0;...
	0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;0;-24;24;-12;0;0;-12;12;-12;0;0;-24;-12;-12;0;48;0;24;0;0;0;-12;0;-12;...
	-12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;24;48;0;12;-12;0;0;-12;0;-12;-12;...
	-12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;24;48;-24;0;0;-12;-12;-12;0;-24;...
	0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48]);
lk(tril(ones(24))==1)=ke';
lk = reshape(lk,24,24 );
lk =lk + lk'- diag( diag( lk ) );
Num_node = (1+nely)*(1+nelx)*(1+nelz);
nodenrs = reshape(1:Num_node,1+nely,1+nelx,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1,nelx*nely*nelz,1);
edofMat = repmat(edofVec,1,24)+ repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);jH = ones(size(iH));sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
	for i1 = 1:nelx
		for j1 = 1:nely
			e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
			for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
				for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
					for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
						e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
						k = k+1;
						iH(k) = e1;
						jH(k) = e2;
						sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
					end
				end
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
	xTilde(1:nely,1:nelx,1:nelz,i+1)=reshape(x(1+i*NDV3:NDV3+i*NDV3),nely,nelx,nelz);
end
xPhys = (tanh(beta*0.5) + tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5) + tanh(beta*(1-0.5)));
%% MMA PARAMETER INITIALIZE
a0=0; a=zeros(NMaterial,1);cc=1.0e6*ones(NMaterial,1); d=0*ones(NMaterial,1);
xmin=1e-3*ones(NDV,1); xmax=1*ones(NDV,1);
xold=x; xolder=xold;low=0; upp=1;
%% START ITERATION
for iter=1:NumIter
	loopbeta=loopbeta+1;
	%% FE-ANALYSIS
	KE0=lk; p=6;
	for i=1:NMaterial
		KE(:,:,i)=par.E(i)*KE0;
	end
	KE_test=zeros(24*24,nely*nelx*nelz); g=zeros(NMaterial,1);
	r1=xPhys(:,:,:,1); r2=xPhys(:,:,:,2); r3=xPhys(:,:,:,3);
	for j=1:NMaterial
		post_x(1:nely,1:nelx,1:nelz,j)=(xPhys(:,:,:,j).*(r1.^p + r2.^p + r3.^p ).^(1/p))./(r1 + r2 + r3+delta);
		dr(1:nely,1:nelx,1:nelz,j)=(post_x(1:nely,1:nelx,1:nelz,j)).^penal;
		KE_test=reshape(KE(:,:,j),24*24,1)*reshape(dr(:,:,:,j),1,nely*nelx*nelz)+KE_test;
		g(j)=sum(sum(sum(xPhys(:,:,:,j))));
	end
	sK = reshape(KE_test,24*24*nelx*nely*nelz,1);
	K = sparse(iK,jK,sK);
	K = (K+K')/2;
	%% SOLVING
	U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
	U(fixeddofs,:)= 0;
	c=F'*U;
	%% COMPUTE SENSITIVITIES
	dc=zeros(nely,nelx,nelz,NMaterial);
	dgdx_test=zeros(NMaterial,nely*nelx*nelz*NMaterial);
	KES=reshape(sum(U(edofMat)*KE0.*U(edofMat),2),[nely,nelx,nelz]);
	for n=1:NMaterial
		for m=1:NMaterial
			if m==n
				drmdn(1:nely,1:nelx,1:nelz,m)=penal.*post_x(1:nely,1:nelx,1:nelz,m).^(penal -1).*...
					((r1.^p + r2.^p + r3.^p).^(1/p)./(r1 + r2 + r3+delta) - (xPhys(:,:,:,m).*(r1.^p + r2.^p + r3.^p).^(1/p))./(r1 + r2 + r3+delta).^2 + (xPhys(:,:,:,m).*xPhys(:,:,:,n).^(p - 1).*(r1.^p + r2.^p + r3.^p).^(1/p - 1))./(r1 + r2 + r3+delta));
				mass_drmdn=1;
			else
				drmdn(1:nely,1:nelx,1:nelz,m)=penal.*post_x(1:nely,1:nelx,1:nelz,m).^(penal - 1).*...
					(((xPhys(:,:,:,m).*xPhys(:,:,:,n).^(p - 1).*(r1.^p + r2.^p + r3.^p).^(1/p - 1))./(r1 + r2 + r3+delta)- (xPhys(:,:,:,m).*(r1.^p + r2.^p + r3.^p).^(1/p))./(r1 + r2 + r3+delta).^2));
				mass_drmdn=0;
			end
			dgdx_test(m,1+(n-1)*NDV3:NDV3+NDV3*(n-1))=repmat(mass_drmdn,1,nelx*nely*nelz);
		end
		dc(:,:,:,n)=-drmdn(:,:,:,1).*KES.*par.E(1)-drmdn(:,:,:,2).*KES.*par.E(2)-drmdn(:,:,:,3).*KES.*par.E(3);
	end
	%% FILTERING/MODIFICATION OF SENSITIVITIES
	filtereddc=[];dfdx=[];
	for i=1:NMaterial
		dxx = beta * (1-tanh(beta*(xTilde-0.5)).*tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5) + tanh(beta*(1-0.5)));
		[filtereddc(:,:,:,i)]= H*(reshape(dc(:,:,:,i).*dxx(:,:,:,i),NDV3,1)./Hs);
		dfdx=[dfdx;filtereddc(:,:,:,i)];
		g(i)=(g(i)/NDV3-0.5/NMaterial)*1;
		dgdx1=(dgdx_test(i,:)/NDV3)*1;
		for j=0:NMaterial-1
			dgdx2=reshape(dgdx1(:,1+j*NDV3:NDV3+j*NDV3),nely,nelx,nelz);
			dgdx(i,1+j*NDV3:NDV3+j*NDV3)=H*(reshape(dgdx2(:,:,:).*dxx(:,:,:,i),NDV3,1)./Hs);
		end
	end
	%% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
	for j=1:length(x)
		xmin(j)=max(1e-3,x(j)-movelimit); xmax(j)=min(1,x(j)+movelimit);
	end
	%% MMA OPTIZATION
	[xnew, y, z, lamda, ksi, eta, mu, zeta, s, low, upp] = ...
		mmasub(NMaterial,length(x), iter, ...
		x, xmin, xmax, xold, xolder, ...
		c, (dfdx), g*1e4, dgdx*1e4, low, upp, a0, a, cc, d);
	%% PRINT RESULTS
	disp(sprintf('Iter.:%3d Obj.: %8.4e max constraint: %6.3e', iter,c ,max(g)))
	%% UPDATE PARAMETER
	change=norm(abs(xnew-x));
	xolder=xold; xold=x; x = xnew;
	for i=0:NMaterial-1
		xTilde(1:nely,1:nelx,1:nelz,i+1)= reshape((H*xnew(1+i*NDV3:NDV3+i*NDV3,:))./Hs,nely,nelx,nelz);
	end
	xPhys = (tanh(beta*0.5) + tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5) + tanh(beta*(1-0.5)));
	if beta < 100 && (loopbeta >= 50)
		beta = 2*beta;loopbeta = 0;
	end
	%% PLOT CURRENT MULTI-MATERIAL DENSITIES
	figure(1);
	subplot(2,2,1); display_3D(post_x,1,NMaterial);axis equal; hold on;title('Material 1')
	subplot(2,2,2); display_3D(post_x,2,NMaterial);axis equal; hold on;title('Material 2')
	subplot(2,2,3); display_3D(post_x,3,NMaterial);axis equal; hold on;title('Material 3')
	subplot(2,2,4); display_3D(post_x,4,NMaterial);axis equal; hold on;title('Material Index')
end
%% PLOT DENSITIES
close all;
display_3D(post_x, NMaterial+1,NMaterial);axis equal; hold on;
%% 3D TOPOLOGY DISPLAY FUNCTION
function display_3D(post_x,NMat,Tolnmaterial)
	if NMat<Tolnmaterial+1
		NM=NMat;Nn=NMat;
	else
		NM=NMat-1;Nn=1;
	end
	for Nm=Nn:NM
		[nely,nelx,nelz] = size(post_x(:,:,:,Nm));
		hx = 1; hy = 1; hz = 1;
		face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
		for k = 1:nelz
			z = (k-1)*hz;
			for i = 1:nelx
				x = (i-1)*hx;
				for j = 1:nely
					y = nely*hy - (j-1)*hy;
					if (post_x(j,i,k,Nm) >0.5)
						vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx;x+hx y-hx z+hx;x+hx y z+hx];
						vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
						if Nm==1
							patch('Faces',face,'Vertices',vert,'FaceColor','r');
						elseif Nm==2
							patch('Faces',face,'Vertices',vert,'FaceColor','b');
						elseif Nm==3
							patch('Faces',face,'Vertices',vert,'FaceColor','g');
						end
					end
				end
			end
		end
	end
	axis equal; axis tight; axis off; box on; view([-45,30]); hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zheng R et al. (2024), Applied Sciences, DOI: 10.3390/app14020657           %
% Title:                                                                      %
%   An Efficient Code for the Multi-Material Topology Optimization            %
%   of 2D/3D Continuum Structures Written in Matlab                           %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

