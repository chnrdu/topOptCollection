%%%%%%%%%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE              Nov 2010 %%%%%%%%%%%%
%%%%%%%%%%%% COMMENTED - OUT BY HAOTIAN_W                    AUGUST 2020 %%%%%%%%%%%%
function top88(nelx,nely,volfrac,penal,rmin,ft)
% ===================================================================================
% 88行程序在99行程序上的主要改进:
%   1) 将for循环语句向量化处理，发挥MATLAB矩阵运算优势;
%   2) 为不断增加数据的数组预分配内存，避免MATLAB花费额外时间寻找更大的连续内存块;
%   3) 尽可能将部分程序从循环体里抽出，避免重复计算;
%   4) 设计变量不再代表单元伪密度，新引入真实密度变量xphys;
%   5) 将原先的所有子程序都集成在主程序里，避免频繁调用;
%   总体上，程序的效率有显著提升（近百倍）、内存占用降低，但是对初学者来说可读性不如99行
% ===================================================================================
% nelx    : 水平方向上的离散单元数;
% nely    : 竖直方向上的离散单元数;
%
% volfrac : 容积率，材料体积与设计域体积之比，对应的工程问题就是"将结构减重到百分之多少";
%
% penal   : 惩罚因子，SIMP方法是在0-1离散模型中引入连续变量x、系数p及中间密度单元，从而将离
%           散型优化问题转换成连续型优化问题，并且令0≤x≤1，p为惩罚因子，通过设定p>1对中间密
%           度单元进行有限度的惩罚，尽量减少中间密度单元数目，使单元密度尽可能趋于0或1;
% 
%           合理选择惩罚因子的取值，可以消除多孔材料，从而得到理想的拓扑优化结果:
%               当penal<=2时     存在大量多孔材料，计算结果没有可制造性;
%               当penal>=3.5时   最终拓扑结果没有大的改变;
%               当penal>=4时     结构总体柔度的变化非常缓慢，迭代步数增加，计算时间延长;
%            
% rmin    : 敏度过滤半径，防止出现棋盘格现象;
 
% ft      : 与99行程序不同的是，本程序提供了两种滤波器，
%           ft=1时进行灵敏度滤波，得到的结果与99行程序一样；ft=2时进行密度滤波。
% ===================================================================================
%% 定义材料属性
% E0杨氏弹性模量;
E0 = 1;
 
% Emin自定义空区域杨氏弹性模量，为了防止出现奇异矩阵，这里和99行程序不同，参见论文公式（1）;
% 不要定义成0，否则就没有意义了;
Emin = 1e-9;
 
% nu泊松比;
nu = 0.3;
 
%% 有限元预处理
% 构造平面四节点矩形单元的单元刚度矩阵KE，详见有限元理论推导
% This is for the plane stress condition,
% that means the stress in the z-direction is zero, which is suitable for thin structures
A11 = [12  3 -6 -3; ...
        3 12  3  0; ...
       -6  3 12 -3; ...
       -3  0 -3 12];
A12 = [-6 -3  0  3; ...
       -3 -6 -3 -6; ...
        0 -3 -6  3; ...
        3 -6  3 -6];
B11 = [-4  3 -2  9; ...
        3 -4 -9  4; ...
       -2 -9 -4 -3; ...
        9  4 -3 -4];
B12 = [ 2 -3  4 -9; ...
       -3  2  9 -2; ...
        4  9  2  3; ...
       -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

% nodenrs存放节点编号，按照列优先的顺序，从1到(1+nelx)*(1+nely);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
 
% edofVec存放所有单元的第一个自由度编号（左下角），参见下面图1;
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
 
% edofMat按照行存放每个单元4个节点8个自由度编号，所以列数是8，行数等于单元个数;
% 存放顺序是：[左下x  左下y  右下x  右下y  右上x  右上y  左上x  左上y];
% 第一个repmat将列向量edofVec复制成8列，其他自由度从第一个自由度上加或减可以得到，参见下面图2;
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
 
% 根据iK、jK和sK三元组生成总体刚度矩阵的稀疏矩阵K，K = sparse(iK,jK,sK)
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
 
% 施加载荷，直接构造稀疏矩阵
% generate 2*(nely+1)*(nelx+1) x 1 sparse matrix F, 
% F is -1 at 2nd DOF of node 1
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
 
% 施加约束，与99行程序相同，唯一的区别是这里把“定义载荷约束”放在了循环体外，提高效率
% Initalize the displacement vector U as a zero vector
U = zeros(2*(nely+1)*(nelx+1),1);
% THE FOLLOWING DOFS (VALUES OF DOF MATRIXES) ARE 
% THE INDEXES OF THE NODES' DOFS IN A VECTOR
% Define the fixed DOFs, which are 
% 1st(x-axis) DOFs of each nodes on the left edge,
% and 2nd (y-axis) DOF of the bottom left corner node
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
% All DOFs
alldofs = [1:2*(nely+1)*(nelx+1)];
% Free DOFs, which are not fixed
freedofs = setdiff(alldofs,fixeddofs);

%% 滤波预处理
% 根据iH、jH和sH三元组生成加权系数矩阵的稀疏矩阵H，H = sparse(iH,jH,sH);
% H是论文公式（8）中的Hei,Hs是论文公式（9）中的Sigma(Hei);
 
% 为数组iH、jH、sH预分配内存;
% improve the efficiency by preallocating memory for arrays iH, jH, sH
total_ele = nelx*nely;
total_neighbors= (2*(ceil(rmin)-1)+1)^2;
total_entries = total_ele * total_neighbors;
iH = ones(total_entries,1);
jH = iH;
sH = zeros(total_entries,1);
e1 = 0;
k = 0;

% 4层for循环，前两层i1、j1是遍历所有单元，后两层i2、j2是遍历当前单元附近的单元
% 这种敏度过滤技术的本质是利用过滤半径范围内各单元敏度的加权平均值代替中心单元的敏度值
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = e1 + 1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));  % 加权系数
      end
    end
  end
end
% 所生成的矩阵H是一个稀疏矩阵，每一行代表一个单元，每一行的元素是该单元与其他单元的加权系数
% e.g., The value of H(i,j) is the weight of the j(e2)-th element in the i(e1)-th element's neighborhood
% where, the i,j is determined by the for-loop above, and the weight(sH) is determined by the distance
%         +---+---+---+............
%         |-x-|-x-|-x-|............
%         +---+---+---+............
%         |-x-|i,j|-x-|............
%         +---+---+---+............
%         |-x-|-x-|-x-|............
%         +---+---+---+............
H = sparse(iH,jH,sH,nelx*nely,nelx*nely);
% Hs is the sum of each row of H
% sum the 2nd dimension of matrix H
Hs = sum(H,2);
%% 迭代初始化
x = repmat(volfrac,nely,nelx);   % x设计变量;
xPhys = x;                       % xphys单元物理密度，真实密度，这里与99行不一样;
loop = 0;                        % loop存放迭代次数;
change = 1; 
 
%% 进入优化迭代，到此为止上面的部分都是在循环外，比99行效率提高很多
while change > 0.001 && loop < 500
  loop = loop + 1;

  %% 有限元分析求解
  % (Emin+xPhys(:)'.^penal*(E0-Emin))就是论文公式（1），由单元密度决定杨氏弹性模量;
  % KE is 8 x 8, E is nelx*nely x 1, total elements in sK is 64*nelx*nely
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);

  % 组装总体刚度矩阵的稀疏矩阵;
  % K = (K+K')/2确保总体刚度矩阵是完全对称阵，因为这会影响到MATLAB求解有限元方程的算法;
  % 当K是实对称正定矩阵时则采用Cholesky平方根分解法，反之则采用速度更慢的LU三角分解法;
  K = sparse(iK,jK,sK); K = (K+K')/2;

  % 正式求解,KU=F;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  
  %% 目标函数和关于真实密度场的灵敏度信息
  % 参见论文公式（2），ce是ue^T*k0*ue，c是目标函数;
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  
  % 参见论文公式（5），只不过将设计变量x换成真实密度xphys;
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  
  % 参见论文公式（6）;
  dv = ones(nely,nelx);
  
  %% 敏度滤波 或 密度滤波
  
  if ft == 1
    % ft=1灵敏度滤波（得到的结果同99行程序），参见论文公式（7）;
    % 完成敏度滤波;
    % 1e-3是公式（7）中的gamma，max(1e-3,x(:))是为了防止分母出现0;
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    
  elseif ft == 2
    % ft=2密度滤波
    % 密度滤波进行两个操作，一个是密度滤波，在下面的循环里面
    % 另一个是根据链式法则修正目标函数和体积约束的灵敏度信息，就是这里，参见公式（10）;
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  
  %% OC优化准则法更新设计变量和单元密度
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    
    % 更新设计变量，参见论文公式（3）;
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    
    if ft == 1
      % 敏度滤波没有密度滤波那么复杂，设计变量就是当前单元的伪密度;
      xPhys = xnew;
      
    elseif ft == 2
      % 完成密度滤波，参见论文公式（9），设计变量经过滤波之后才是单元伪密度;
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid;
    else 
      l2 = lmid; 
    end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;

  %% 显示结果（同99行程序）
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  colormap(gray); imagesc(1-xPhys); clim([0 1]); axis equal; axis off; pause(1e-3);
end
%————————————————
%
%版权声明：本文为博主原创文章，遵循 CC 4.0 BY-SA 版权协议，转载请附上原文出处链接和本声明。
%                        
%原文链接：https://blog.csdn.net/BAR_WORKSHOP/article/details/108287668
% other references:
% https://blog.qiql.net/archives/fem1
% https://abg.baidu.com/view/39bacd212f60ddccda38a0eb
% https://www.jishulink.com/post/1885570
% https://enterfea.com/fea-by-hand-2d-plates/
% https://engcourses-uofa.ca/books/introduction-to-solid-mechanics/finite-element-analysis/two-dimensional-solid-elements/quadrilateral-elements/

% Efficient topology optimization in MATLAB using 88 lines of code
% by Erik Andreassen, Anders Clausen, Mattias Schevenels, 
%    Boyan S. Lazarov, Ole Sigmund
% Struct Multidisc Optim (2011) 43:1–16
% DOI 10.1007/s00158-010-0594-7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%