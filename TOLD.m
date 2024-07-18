function TOLD
%initialize---system parameters
dbstop if error
loop=0;wd=cd;
%initialize---material parameters
E0=1e8;nu=0.3;penal=3;c2=5e-4;strain0=0.6;
%initialize---boundary condition and mesh
[nele,nelx,nely,thickness,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter;
%initialize---filiter
rmin=1.8;
[H,Hs]=preparefilter(rmin,nele,nelx,nely);
volfrac=0.25;
x = repmat(volfrac,[nele,1]);
xba = (H*x(:))./Hs;  
%iinitialize--MMA parameters
changeobj=1;m=1;a0=1;amma=0;cmma=10000;
d=0;xold1=x;xold2=x;dcdxba(nele,1)=0;
xmax=ones(nele,1);xmin  = 0.001*ones(nele,1);
low=xmin;upp=xmax;
tt=zeros(3,100);
while changeobj>1e-4 && loop < 100   
    tic
    loop = loop+1;   
    %generatecommand 
    generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,edofMat,fixeddof,loaddof,penal,c2,Kin,Kout,outdof);   
    tt(1,loop)=toc;
    %run ANSYS APDL as a subroutine
    cd 'D:\software\ansys14\ANSYS Inc\v140\ansys\bin\winx64'
    !ANSYS140 -b -p ane3fl -i D:\command.txt -o D:\1.out  
    tt(2,loop)=toc;
    cd(wd)    
    %read nodal displacement solution
    U=load('D:\NODALDISPLACEMENT.txt');
    delete('D:\NODALDISPLACEMENT.txt');
    U=U';U=U(:);       
    obj(loop)=U(outdof);
    %read the element logarithmic equivilant strain and update c2
    strain=load('D:\eqvstrain.txt');
    strain=reshape(strain,2,[]);
    maxvms=max(strain');svon=maxvms(2);
    if svon<=strain0
        c2=max((svon/strain0)^0.5,0.8)*c2;
    else
        c2=(svon/strain0)^2*c2;        
    end
    c2=max(c2,1e-6);
    %read the element strain energy 
    strain_energy=load('D:\strainenergy.txt');
    strain_energy99=load('D:\strainenergy99.txt');    
    %sensitivity analysis
    dstrain_energy=reshape(strain_energy99-strain_energy,2,[]);    
    for i=1:nele
    dcdxba(i)=-1/(F(loaddof)/100)*(dstrain_energy(1,i)*penal/xba(i)-...
        dstrain_energy(2,i)*penal*xba(i)^(penal-1)/(1-xba(i)^penal+1e-10));
    end
    dvdxba=ones(nely,nelx);         
    dcdx = H*(dcdxba(:)./Hs);                                                     
    dvdx = H*(dvdxba(:)./Hs);      
    %save data
    lp=num2str(loop);       
    save(lp);         
    %update the design variables by MMA
    ratemma=11000;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
    mmasub(m,nele,loop,x,xmin,xmax,xold1,xold2,ratemma(1)*obj(loop),ratemma(1)*dcdx,...
    0*dcdx,sum(xba)-nele*volfrac,dvdx',0*dvdx',low,upp,a0,amma,cmma,d);    
    xold2 = xold1;
    xold1 = x;
    x = xmma;  
    xba(:) = (H*x(:))./Hs;       
    if loop>=6
    changeobj=mean((abs(obj(loop-4:loop)-obj(loop-5:loop-1))))/abs(obj(loop));
    end
    fprintf('It.:%5i Obj.:%11.4f Vol.:%2.4f ch.:%2.5f\n',loop,full(obj(loop)),mean(xba),changeobj);          
    colormap(gray); imagesc(1-reshape(xba,nely,nelx)); axis equal; axis tight; axis off;pause(1e-6);  
    tt(3,loop)=toc;
end
plot(1:loop,obj);
end
function [nele,nelx,nely,cc,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter
nelx=100;nely=50;
aa=100e-3/nelx;bb=50e-3/nely;cc=1e-3;
[i0,j0] = meshgrid(0:nelx,0:nely);
gobalcoords0x=i0*aa;
gobalcoords0y=bb*nely-j0*bb;
gobalcoords0xy=[gobalcoords0x(:) gobalcoords0y(:)]';
gobalcoords0=gobalcoords0xy(:);
[il,jl] = meshgrid(0,nely);
loadnid = il*(nely+1)+(nely+1-jl); 
loaddof = 2*loadnid(:)-1; 
force=5;
Kin=500;
[ilo,jlo] = meshgrid(nelx,nely);
outnid = ilo*(nely+1)+(nely+1-jlo);  
outdof = 2*outnid(:)-1; 
Kout=100;
[iif,jf] = meshgrid(0,0:nely/5);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); 2*fixednid(:)-1];
[iif,jf] = meshgrid(0:nelx,nely);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); fixeddof];
nele = nelx*nely;
ndof = 2*(nelx+1)*(nely+1);
F = sparse(loaddof,1,force,ndof,1);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);                               
edofVec = 2*nodeids(:)+1;                                                              
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2*(nely+1)-2 2*(nely+1)-1 -2 -1],nele,1);%%%%b=[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],����ڵ���ȡ��д������������Ե�Ԫ��1�ڵ�x�������ɶȵ�λ��
end
function [H,Hs]=preparefilter(rmin,nele,nelx,nely)
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
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
end
function  generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,edofMat,fixeddof,loaddof,penal,c2rate,Kin,Kout,outdof);    
fid = fopen('D:\command.txt','w');
ndot=ndof/2;
fprintf(fid,'/CWD,''D:\\''\n');
fprintf(fid,'/PREP7\n');
%define node
for i=1:ndot
fprintf(fid,'N,%d,%G,%G,%G\n',i,gobalcoords0(2*i-1),gobalcoords0(2*i),0);
end
%define element type and the thickness
fprintf(fid,'et,1,plane182\nKEYOPT,1,1,0\nKEYOPT,1,3,3\nkeyopt,1,6,0\nR,1,%G, \nMPTEMP,1,0 \n',thickness);
edotMat=0.5*edofMat(:,[2,4,6,8]);
%calculate the element parameters
MU0=E0*xba.^penal/(2*(1+nu));
K0=E0*xba.^penal/(2*(1-2*nu));
d=2./K0;
c1=(1-xba.^penal)*E0*1e-9/6;
c2=(1-xba.^penal)*E0*c2rate;
for i=1:nele
        %define original elements
        fprintf(fid,'TB,HYPE,%d,1,2,NEO\nTBTEMP,0\n',2*i-1);
        fprintf(fid,'TBDATA,,%G,%G,,,,\n',MU0(i),d(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i-1);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
        %define the additive hyperelastic elements  
        fprintf(fid,'TB,HYPE,%d,1,2,YEOH\nTBTEMP,0\n',2*i);
        fprintf(fid,'TBDATA,,%G,%G,1e-9,1e6,,\n',c1(i),c2(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
end
%%%%%define node of the spring
coords0w=reshape(gobalcoords0,2,[]);
lx=max(coords0w(1,:))-min(coords0w(1,:));
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+1,gobalcoords0(loaddof)+2*lx,gobalcoords0(loaddof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+1);
fprintf(fid,'d,%d,uy,0\n',ndot+1);
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+2,gobalcoords0(outdof)+2*lx,gobalcoords0(outdof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+2);
fprintf(fid,'d,%d,uy,0\n',ndot+2);
%%%%%define the spring
fprintf(fid,'ET,2,LINK180\nKEYOPT,2,2,0\nR,2,1, ,0\n');
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+1,Kin*2*lx,2*nele+1,0.3,2*nele+1,(loaddof+1)/2,ndot+1);
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+2,Kout*2*lx,2*nele+2,0.3,2*nele+2,(outdof+1)/2,ndot+2);
%apply the displacement
nfix=size(fixeddof,1);
for i=1:nfix
    if mod(fixeddof(i),2)==1
        fprintf(fid,'d,%d,ux,0\n',(fixeddof(i)+1)/2);
    else       
        fprintf(fid,'d,%d,uy,0\n',fixeddof(i)/2);
    end
end
%apply the external load
nload=size(loaddof,1);
for i=1:nload
    if mod(loaddof(i),2)==1
    fprintf(fid,'F,%d,fx,%G\n',(loaddof(i)+1)/2,full(F(loaddof(i))));
    else
    fprintf(fid,'F,%d,fy,%G\n',loaddof(i)/2,full(F(loaddof(i))));
    end
end
%solve 
fprintf(fid,'finish\n/sol\nANTYPE,0\nNLGEOM,1\nNSUBST,1,0,0\n');
fprintf(fid,'CNVTOL,U,-1, \n'); 
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',sum(abs(full(F))));
fprintf(fid,'OUTRES,ERASE\nOUTRES,ALL,ALL\n/status,solu\nsolve\n');
%apply the addictional external force
fprintf(fid,'F,%d,fx,%f\n',(outdof+1)/2,full(F(loaddof)/100));
fprintf(fid,'NSUBST,1,0,0\n');
%��������׼��
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',1e4*sum(abs(full(F))));
fprintf(fid,'solve\nfinish\n');
%post processing---write the element strain energy 
fprintf(fid,'/POST1\n');
fprintf(fid,'SET, , ,1, ,1, , \n');
fprintf(fid,'esel,all\n');
fprintf(fid,'*dim,STEN,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea,SENE\n');
fprintf(fid,'*vget,STEN,elem,1,etab,ea\n');
fprintf(fid,'*cfopen,strainenergy,txt\n');
fprintf(fid,'*vwrite,STEN(1,1)\n');
fprintf(fid,'(E13.5)\n');
%post processing---the element logarithmic equivilant strain
fprintf(fid,'*dim,strain1,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ab,EPEL,EQV\n');
fprintf(fid,'*vget,strain1,elem,1,etab,ab\n');
fprintf(fid,'*cfopen,eqvstrain,txt\n');
fprintf(fid,'*vwrite,strain1(1,1)\n');
fprintf(fid,'(E13.5)\n');
%post processing---the nodal displacement
fprintf(fid,'nsel,all\n');
fprintf(fid,'*dim,UARRAY,array,%d,2\n',ndot);
fprintf(fid,'*vget,UARRAY(1,1),node,1,u,x\n');
fprintf(fid,'*vget,UARRAY(1,2),node,1,u,y\n');
fprintf(fid,'*cfopen,NODALDISPLACEMENT,txt\n');
fprintf(fid,'*vwrite,UARRAY(1,1),UARRAY(1,2)\n');
fprintf(fid,'(E13.5,5X,E13.5)\n');
%post processing---write the element strain energy 
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