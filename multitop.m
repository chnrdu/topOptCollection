function multitop(nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf)
%%  115 LINES MATLAB CODE MULTIPHASE MINIMUM COMPLIANCE TOPOLOGY OPTIMIZATION
%%  nx,ny: number of elements in x and y directions
%%  tol_out: tolerance for outer iterations
%%  tol_f: tolerance for filter update
%%  iter_max_in: maximum number of inner iterations
%%  iter_max_out: maximum number of outer iterations
%%  p: number of phases
%%  q: penalty factor
%%  e: Young's modulus vector of phases
%%  v: volume fraction vector of phases
%%  rf: filter radius
  alpha = zeros(nx*ny,p);     % initialize design variables with domain size
  for i = 1:p                 % set variables with different phases
    alpha(:,i) = v(i);
  end
  %% MAKE FILTER
  [H,Hs] = make_filter(nx,ny,rf);
  change_out = 2*tol_out; iter_out = 0;
  while (iter_out < iter_max_out) && (change_out > tol_out)
    alpha_old = alpha;
    for a = 1:p
      for b = a+1:p
        [obj,alpha] = bi_top(a,b,nx,ny,p,v,e,q,alpha,H,Hs,iter_max_in);
      end
    end
    iter_out = iter_out + 1;
    change_out = norm(alpha(:)-alpha_old(:),inf);
    fprintf('Iter:%5i Obj.:%11.4f change:%10.8f\n',iter_out,obj,change_out);
    %% UPDATE FILTER
    if (change_out < tol_f) && (rf>3)
      tol_f = 0.99*tol_f; rf = 0.99*rf; [H,Hs] = make_filter (nx,ny,rf);
    end
    %% SCREEN OUT TEMPORAL TOPOLOGY EVERY 5 ITERATIONS
    if mod(iter_out,5)==0
      I = make_bitmap (p,nx,ny,alpha);
      image([flipdim(I ,2) I]), axis image off, drawnow;
    end
  end
end
%% MAKE FILTER
function [H,Hs] = make_filter (nx,ny,rmin)
  ir = ceil(rmin)-1;
  iH = ones(nx*ny*(2*ir+1)^2,1);
  jH = ones(size(iH)); sH = zeros(size(iH)); k = 0;
  for i1 = 1:nx
    for j1 = 1:ny
      e1 = (i1-1)*ny+j1;
      for i2 = max(i1-ir,1):min(i1+ir,nx)
        for j2 = max(j1-ir,1):min(j1+ir,ny)
          e2 = (i2-1)*ny+j2; k = k+1; iH(k) = e1; jH(k) = e2;
          sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
        end
      end
    end
  end
  H = sparse(iH,jH,sH); Hs = sum(H,2);
end
%% MODIFIED BINARY-PHASE TOPOLOGY OPTIMIZATION SOLVER
function [o,alpha] = bi_top(a,b,nx,ny,p,v,e,q,alpha_old,H,Hs,iter_max_in)
  alpha = alpha_old; iter_in = 0; nu = 0.3;
  %% PREPARE FINITE ELEMENT ANALYSIS
  % 2D plane stress constitutive matrix
  A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
  A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
  B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
  B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
  KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
  % element degrees of freedom, node matrix indices
  % reshape 1:(1+nx)*(1+ny) into 1+ny by 1+nx matrix
  nodenrs = reshape(1:(1+nx)*(1+ny),1+ny,1+nx);
  % The vector for the 1st element DOF (x-axis)(lower left corner of element),
  % e.g., 3 is the value for element 1, shown in the figure below
  % ny by nx matrix to 1D vector
  edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nx*ny,1);
  % Generate the element DOF matrix for all elements
  % 8 DOs per element (2 DOF per node, 4 nodes per element)
  % This is the DOF index for an element,
  % e.g., edofMat(1,:) is the DOF index for element 1
  %
  % 1,2 ..., (2*ny+3),2*ny+4 
  % |                |
  % 3,4 ..., (2*ny+5),2*ny+6
  edofMat = repmat(edofVec,1,8)+repmat([0 1 2*ny+[2 3 0 1] -2 -1],nx*ny,1);
  % outer loop
  % nx*ny*8 by 8*1 sparse matrix to 64*nx*ny by 1 vector
  iK = reshape(kron(edofMat,ones(8,1))',64*nx*ny,1);
  jK = reshape(kron(edofMat,ones(1,8))',64*nx*ny,1);
  %% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
  % specify F as 2*(ny+1)(nyx+1) by 1 sparse matrix
  % with F(2,1) = -1
  F = sparse(2,1,-1,2*(ny+1)*(nx+1),1); % MBB
  % for MBB, set the indices of x- of central cutting line and 
  % y- of the lower right corner as fixed
  fixeddofs = union([1:2:2*(ny+1)],[2*(nx+1)*(ny+1)]);
  U = zeros(2*(ny+1)*(nx+1),1);
  alldofs = [1:2*(ny+1)*(nx+1)];
  freedofs = setdiff(alldofs,fixeddofs);
  %% INNER ITERATIONS
  while iter_in < iter_max_in
    iter_in = iter_in + 1;
    %% FE-ANALYSIS
    E = e(1)*alpha(:,1).^q;                              % Can be changed %
    for phase = 2:p
      E = E + e(phase)*alpha(:,phase).^q;                % Can be changed %
    end
    size(E)
    sK = reshape(KE(:)*E(:)',64*nx*ny,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = sum((U(edofMat)*KE).*U(edofMat),2);             % Can be changed %
    o = sum(sum(E.*ce));
    dc = -(q*(e(a)-e(b))*alpha(:,a).^(q-1)).*ce;         % Can be changed %
    %% FILTERING OF SENSITIVITIES
    dc = H*(alpha(:,a).*dc)./Hs./max(1e-3,alpha(:,a)); dc = min(dc,0);
    %% UPDATE LOWER AND UPPER BOUNDS OF DESIGN VARIABLES
    move = 0.2;
    r = ones(nx*ny,1);
    for k = 1:p
      if (k ~= a) && (k ~= b)
        r = r - alpha(:,k);
      end
    end
    l = max(0,alpha(:,a)-move);
    u = min(r,alpha(:,a)+move);
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
    l1 = 0; l2 = 1e9;
    while (l2-l1)/(l1+l2) > 1e-3
      lmid = 0.5*(l2+l1);
      alpha_a = max(l,min(u,alpha(:,a).*sqrt(-dc./lmid)));
      if sum(alpha_a) > nx*ny*v(a); l1 = lmid; else l2 = lmid; end
    end
    alpha(:,a) = alpha_a;
    alpha(:,b) = r-alpha_a;
  end
end
%% MAKE BITMAP IMAGE OF MULTIPHASE TOPOLOGY
function I = make_bitmap (p,nx,ny,alpha)
  color = [1 0 0; 0 0 .45; 0 1 0; 0 0 0; 1 1 1];
  I = zeros(nx*ny,3);
  for j = 1:p
    I(:,1:3) = I(:,1:3) + alpha(:,j)*color(j,1:3);
  end
  I = imresize(reshape(I,ny,nx,3),10,'bilinear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is from the paper                                               %
% "Alternating active-phase algorithm for multimaterial                    %
%  topology optimization problems: a 115-line MATLAB implementation"       %
% by Rouhollah Tavakoli & Seyyed Mohammad Mohseni  (2014),                 %
% Structural and Multidisciplinary Optimization,                           %
% Vol 49, pp.621--642.  DOI:10.1007/s00158-013-0999-1                      %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
