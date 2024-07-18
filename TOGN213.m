%A 213-line topology optimization code for geometrically nonlinear structures
function TOGN213
%1 Initialization
dbstop if error
loop=0;wd=cd;
E0=3e9;nu=0.4;penal=3;c2=5e-4;strain0=0.6;
[nele,nelx,nely,thickness,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof]=initialize2Dcantilever;
rmin=nely/10;
[H,Hs]=preparefilter(rmin,nele,nelx,nely);
volfrac=0.8;
x = repmat(volfrac,[nele,1]);
xba = (H*x(:))./Hs;  
changeobj=1;m=1;a0=1;amma=0;cmma=10000;
d=0;xold1=x;xold2=x;dcdxba(nele,1)=0;
xmax=ones(nele,1);xmin=0.001*ones(nele,1);
low=xmin;upp=xmax;
%start iteration
while changeobj>1e-4 && loop < 100      
    loop = loop+1;
    volfrac=max(volfrac-0.3/2,0.5);
    %2 Commands generation
    generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,edofMat,fixeddof,loaddof,penal,c2);    
    %3 Design variable updating    
    %3.1 run ANSYS APDL as a subroutine
    cd 'D:\software\ansys14\ANSYS Inc\v140\ansys\bin\winx64'
    !ANSYS140 -b -p ane3fl -i D:\command.txt -o D:\1.out  
    cd(wd)        
    %3.2 read results    
    U=load('D:\NODALDISPLACEMENT.txt');    
    delete('D:\NODALDISPLACEMENT.txt');
    U=U';U=U(:);       
    obj(loop)=F'*U;
    strain=load('D:\eqvstrain.txt');
    strain=reshape(strain,2,[]);
    maxvms=max(strain');svon=maxvms(2);
    if svon<=strain0
        c2=max((svon/strain0)^0.5,0.8)*c2;
    else
        c2=(svon/strain0)^2*c2;        
    end
    c2=max(c2,1e-6);
    strain_energy=load('D:\strainenergy.txt');
    strain_energy99=load('D:\strainenergy99.txt');    
    %3.3 sensitivity analysis
    dstrain_energy=reshape(strain_energy99-strain_energy,2,[]);    
    for i=1:nele
        dcdxba(i)=100*(dstrain_energy(1,i)*penal/xba(i)-...
        dstrain_energy(2,i)*penal*xba(i)^(penal-1)/(1-xba(i)^penal+1e-10));
    end
    dvdxba=ones(nely,nelx);         
    dcdx = H*(dcdxba(:)./Hs);                                                     
    dvdx = H*(dvdxba(:)./Hs);          
    lp=num2str(loop);save(lp);         
    %3.4 update the design variables by MMA
    ratemma=10/abs(obj(1));
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
    mmasub(m,nele,loop,x,xmin,xmax,xold1,xold2,ratemma(1)*obj(loop),ratemma(1)*dcdx,...
    0*dcdx,sum(xba)-nele*volfrac,dvdx',0*dvdx',low,upp,a0,amma,cmma,d);    
    xold2 = xold1;xold1 = x;x = xmma;  
    xba(:) = (H*x(:))./Hs;       
    if loop>=6
        changeobj=mean((abs(obj(loop-4:loop)-obj(loop-5:loop-1))))/obj(loop);
    end
    %3.5 plot results
    fprintf('It.:%5i Obj.:%11.4f Vol.:%2.4f ch.:%2.5f\n',loop,full(obj(loop)),mean(xba),changeobj);          
    colormap(gray); imagesc(1-reshape(xba,nely,nelx)); axis equal; axis tight; axis off;pause(1e-6);  
end
end
function [nele,nelx,nely,cc,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof]=initialize2Dcantilever
%mesh 
nelx=120;nely=30;
aa=0.12/nelx;bb=0.03/nely;cc=0.001;
nele = nelx*nely;
ndof = 2*(nelx+1)*(nely+1);
%calculate the coordinate of nodes
[i0,j0] = meshgrid(0:nelx,0:nely);
gobalcoords0x=i0*aa;
gobalcoords0y=bb*nely-j0*bb;
gobalcoords0xy=[gobalcoords0x(:) gobalcoords0y(:)]';
gobalcoords0=gobalcoords0xy(:);
%initialize the boundary condition
[il,jl] = meshgrid(nelx-2:nelx,nely/2-1:nely/2+1);
loadnid = il*(nely+1)+(nely+1-jl); 
loaddof = 2*loadnid(:); 
force=-500/9*ones(9,1);
F = sparse(loaddof,1,force,ndof,1);
[iif,jf] = meshgrid(0,0:nely);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); 2*fixednid(:)-1];
%get the matrix :edofMat
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1); 
edofVec = 2*nodeids(:)+1;
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2*(nely+1)-2 2*(nely+1)-1 -2 -1],nele,1);
end
function [H,Hs]=preparefilter(rmin,nele,nelx,nely)
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1); jH = ones(size(iH)); sH = zeros(size(iH));k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
          for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
               for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                   e2 = (i2-1)*nely+j2;
                   k = k+1; iH(k) = e1; jH(k) = e2;
                   sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                end
           end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
end
function generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,edofMat,fixeddof,loaddof,penal,c2rate)
fid = fopen('D:\command.txt','w');
ndot=ndof/2;
fprintf(fid,'/CWD,''D:\\''\n');
fprintf(fid,'/PREP7\n');
%2.1 modeling--2.1.1 define node
for i=1:ndot
fprintf(fid,'N,%d,%G,%G,%G\n',i,gobalcoords0(2*i-1),gobalcoords0(2*i),0);
end
%2.1.2 define element--difine element type and the thickness
fprintf(fid,'et,1,plane182\nKEYOPT,1,1,0\nKEYOPT,1,3,3\nkeyopt,1,6,0\nR,1,%G, \nMPTEMP,1,0 \n',thickness);
edotMat=0.5*edofMat(:,[2,4,6,8]);
%2.1.2 define element--calculate the equivalent element parameters
MU0=E0*xba.^penal/(2*(1+nu));
K0=E0*xba.^penal/(2*(1-2*nu));
d=2./K0;
c1=(1-xba.^penal)*E0*1e-9/6;
c2=(1-xba.^penal)*E0*c2rate;
for i=1:nele
        %2.1.2 define element--define original elements
        fprintf(fid,'TB,HYPE,%d,1,2,NEO\nTBTEMP,0\n',2*i-1);
        fprintf(fid,'TBDATA,,%G,%G,,,,\n',MU0(i),d(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i-1);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
        %2.1.2 define element--define the additive hyperelastic elements  
        fprintf(fid,'TB,HYPE,%d,1,2,YEOH\nTBTEMP,0\n',2*i);
        fprintf(fid,'TBDATA,,%G,%G,1E-9,1E6,,\n',c1(i),c2(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
end
%2.2.1 boundary condiction of the first load case--apply the displacement
nfix=size(fixeddof,1);
for i=1:nfix
    if mod(fixeddof(i),2)==1
        fprintf(fid,'d,%d,ux,0\n',(fixeddof(i)+1)/2);
    else       
        fprintf(fid,'d,%d,uy,0\n',fixeddof(i)/2);
    end
end
%2.2.2 boundary condiction of the first load case --apply the external load
nload=size(loaddof,1);
for i=1:nload
    if mod(loaddof(i),2)==1
    fprintf(fid,'F,%d,fx,%G\n',(loaddof(i)+1)/2,full(F(loaddof(i))));
    else
    fprintf(fid,'F,%d,fy,%G\n',loaddof(i)/2,full(F(loaddof(i))));
    end
end
%2.2.3 solve the first load case
fprintf(fid,'finish\n/sol\nANTYPE,0\nNLGEOM,1\nNSUBST,5,0,0\n');
fprintf(fid,'CNVTOL,U,-1,\n'); 
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',sum(abs(full(F))));
fprintf(fid,'OUTRES,ERASE\nOUTRES,ALL,ALL\n/status,solu\nsolve\n');
%2.3.1 boundary condiction of the second load case --apply the external load
for i=1:nload
    if mod(loaddof(i),2)==1
    fprintf(fid,'F,%d,fx,%G\n',(loaddof(i)+1)/2,0.99*full(F(loaddof(i))));
    else
    fprintf(fid,'F,%d,fy,%G\n',loaddof(i)/2,0.99*full(F(loaddof(i))));
    end
end
%2.3.2 solve the second load case
fprintf(fid,'NSUBST,1,0,0\n');
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',sum(abs(full(F))));
fprintf(fid,'solve\nfinish\n');
%2.4 post processing--2.4.1 write the element strain energy 
fprintf(fid,'/POST1\n');
fprintf(fid,'SET, , ,1, ,1, , \n');
fprintf(fid,'esel,all\n');
fprintf(fid,'*dim,STEN,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea,SENE\n');
fprintf(fid,'*vget,STEN,elem,1,etab,ea\n');
fprintf(fid,'*cfopen,strainenergy,txt\n');
fprintf(fid,'*vwrite,STEN(1,1)\n');
fprintf(fid,'(E13.5)\n');
%2.4.2 write the element logarithmic equivalent strain
fprintf(fid,'*dim,strain1,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ab,EPEL,EQV\n');
fprintf(fid,'*vget,strain1,elem,1,etab,ab\n');
fprintf(fid,'*cfopen,eqvstrain,txt\n');
fprintf(fid,'*vwrite,strain1(1,1)\n');
fprintf(fid,'(E13.5)\n');
%2.4.3 write  the nodal displacement
fprintf(fid,'nsel,all\n');
fprintf(fid,'*dim,UARRAY,array,%d,2\n',ndot);
fprintf(fid,'*vget,UARRAY(1,1),node,1,u,x\n');
fprintf(fid,'*vget,UARRAY(1,2),node,1,u,y\n');
fprintf(fid,'*cfopen,NODALDISPLACEMENT,txt\n');
fprintf(fid,'*vwrite,UARRAY(1,1),UARRAY(1,2)\n');
fprintf(fid,'(E13.5,5X,E13.5)\n');
%2.4.4 write the element strain energy in the second load case
fprintf(fid,'SET, , ,1, ,2, , \n');
fprintf(fid,'*dim,STEN99,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea99,SENE\n');
fprintf(fid,'*vget,STEN99,elem,1,etab,ea99\n');
fprintf(fid,'*cfopen,strainenergy99,txt\n');
fprintf(fid,'*vwrite,STEN99(1,1)\n');
fprintf(fid,'(E13.5)\n');
fprintf(fid,'finish\n');
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Qi Chen, Xianmin Zhang and Benliang Zhu      %
% Guangdong Key Laboratory of Precision Equipment and Manufacturing Technology,%
% South China University of Technology                                         %
% Please sent your comments to: chenqiscutedu@qq.com                           %
% The code is intended for educational purposes and theoretical details are    %
% discussed in the paper                                                       %
% A 213-line topology optimization code for geometrically nonlinear structures %
% Qi Chen, Xianmin Zhang and Benliang Zhu,                                     %
% Struct Multidisc Optim, 2018                                                 %
% Disclaimer:                                                                  %
% The authors reserves all rights but do not guaranty that the code is free    %
% from errors. Furthermore, we shall not be liable in any event caused by the  %
% use of the program.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%