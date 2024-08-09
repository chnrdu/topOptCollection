%%%% A SOFT-KILL BESO CODE BY X. HUANG and Y.M. XIE %%%%
function sbeso(nelx,nely,volfrac,er,rmin)
    % INITIALIZE
    x(1:nely,1:nelx) = 1.; vol=1.; i = 0; change = 1.; penal = 3.;
    % START iTH ITERATION
    while change > 0.001
        i = i + 1; vol = max(vol*(1-er),volfrac);
        if i > 1; olddc = dc; end
        % FE-ANALYSIS
        [U]=FE(nelx,nely,x,penal);
        % OBJECTIVE FUNCTION AND SENSITTVITY ANALYSIS
        KE=lk;
        c(i)=0.;
        for ely = 1:nely
            for elx = 1:nelx
                n1 = (nely+1) * (elx- 1) + ely; 
                n2 = (nely+1) *  elx     + ely;
                Ue =U([2*n1-1;2*n1; 2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2],1); 
                c(i) = c(i) + 0.5 * x(ely,elx)^penal*Ue'*KE*Ue; 
                dc(ely,elx) = 0.5 * x(ely,elx)^(penal-1)*Ue'*KE*Ue;
            end
        end
        % FILTERING OF SENSITIVITIES
        [dc] =check(nelx,nely,rmin,x,dc);
        % STABLIZATION OF EVOLUTIONARY PROCESS
        if i > 1; dc =(dc+olddc)/2.; end
        % BESO DESIGN UPDATE
        [x] = adddel(nelx,nely,vol,dc,x);
        % PRINT RESULTS
        if i > 10
            change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
        end
        disp([' It.: ' sprintf('%4i',i) ' Obj.: ' sprintf('%10.4f',c(i)) ... 
              'Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely))... 
              ' ch.: ' sprintf('%6.3f',change)])
        % PLOT DENSITIES
        colormap(gray); imagesc(-x); axis equal; axis tight; axis off; pause(1e-6);
    end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=adddel(nelx,nely,volfra,dc,x)
    lim1 = min(min(dc)); lim2 = max(max(dc));
    while ((lim2-lim1)/lim2>1.0e-5)
        th = (lim1+lim2)/2.0;
        x = max(0.001,sign(dc-th));
        if sum(sum(x))-volfra*(nelx*nely) > 0
            lim1=th;
        else
            lim2=th;
        end
    end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcf]=check(nelx,nely,rmin,x,dc)
    dcf=zeros(nely,nelx);
    for i=1:nelx
        for j = 1:nely 
            sum=0.0;
            for m = max(i-floor(rmin),1):min(i+floor(rmin),nelx) 
                for n = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                    fac = rmin - sqrt((i-m)^2 + (j-n)^2);
                    sum = sum+max(0, fac);
                    dcf (j,i) = dcf(j,i) + max(0,fac) * dc(n,m);
                end
            end
            dcf(j,i) =  dcf(j,i)/sum;
        end
    end
%%%%%%%%%% FE- ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx, nely, x, penal)
    [KE] = lk;
    K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
    F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1) ;
    for elx = 1:nelx
        for ely = 1:nely
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)* elx    + ely;
            edof= [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
            K(edof, edof)= K(edof, edof) + x(ely, elx)^penal*KE;
        end
    end
    % DEFINE LOADS AND SUPPORTS (Cantilever)
    F(2*(nelx+1)*(nely+1)-nely, 1) =-1.0;
    fixeddofs= [1:2*(nely+1)] ;
    alldofs = [1:2*(nely+1)*(nelx+1)];
    freedofs = setdiff(alldofs, fixeddofs);
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    U(fixeddofs,:) =  0;
%%%%%%%%%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
    E=1.0;
    nu = 0.3;
    k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
    -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
    KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  This software is published with the following article:
%%  Xiaodong Huang and Yi-Min Xie. A further review of ESO type methods for 
%% topology optimization, Structural and Multidisciplinary Optimization, 
%% 2010, 41(5): 671-683
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
