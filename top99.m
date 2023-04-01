%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Originally commented by Xiaoyu ZHUANG, 2022-11-15
%%                           https://www.zhihu.com/question/61259261
%%   Ref: HAOTIAN_W, TopOpt | 99行拓扑优化程序完全注释. 
%%        https://blog.csdn.net/BAR_WORKSHOP/article/details/108274360
%%        retrieved on 2022-12-04
%%
%%   Comments updated by Run DU <rdu@swjtu.edu.cn>
%%
%%   nelx - x 方向上网格的数目
%%   nely - y 方向上网格的数目
%%   (nely,nelx) matrix
%%   The symbol "+" represents a node.
%%
%%   +-----+-----+---------+-----+-----+-----+-> +---> x
%%   |(1,1)|     |         |     |     |(1,x)|   |
%%   +-----+-----+---------+-----+-----+-----+   |
%%   |(2,1)|     |         |     |     |(2,x)|   V  
%%   +-----+-----+---------+-----+-----+-----+   y
%%   |(3,1)|     |         |     |     |(3,x)|     
%%   +-----+-----n1--------n2----+-----+-----+  n1=(nely+1)*(elx-1)+ely
%%   |     |     |(ely,elx)|     |     |     |  n2=(nely+1)* elx   +ely
%%   +-----+-----+---------+-----+-----+-----+
%%   |     |     |         |     |     |     |     
%%   +-----+-----+---------+-----+-----+-----+
%%   |(y,1)|     |         |     |     |(y,x)|
%%   +-----+-----+---------+-----+-----+-----+
%%   |
%%   V
%%
%%   ^
%%   | Negative load in y-direction
%%
%% ===================================================================================
%% volfrac : 容积率，材料体积与设计域体积之比，对应的工程问题就是"将结构减重到百分之多少";
%%
%% penal   : 惩罚因子，SIMP 方法是在 0-1 离散模型中引入连续变量 x、系数 p 及中间密度单元，从而将离
%%           散型优化问题转换成连续型优化问题，并且令 0≤x≤1，p 为惩罚因子，通过设定 p>1 对中间密
%%           度单元进行有限度的惩罚，尽量减少中间密度单元数目，使单元密度尽可能趋于 0 或 1;
%% 
%%           合理选择惩罚因子的取值，可以消除多孔材料，从而得到理想的拓扑优化结果:
%%               当 penal<=2   时,     存在大量多孔材料，计算结果没有可制造性;
%%               当 penal>=3.5 时,     最终拓扑结果没有大的改变;
%%               当 penal>=4   时,     结构总体柔度的变化非常缓慢，迭代步数增加，计算时间延长;
%%            
%% rmin    : 敏度过滤半径，防止出现棋盘格现象;
%% ===================================================================================
%% 结构优化的数学模型常用名词:
%%   1) 设计变量   在设计中可调整的、变化的基本参数，本算例是单元密度;
%%   2) 目标函数   设计变量的函数，优化设计的目标，本算例是柔度最小;
%%   3) 约束条件   几何、应力、位移等，难点在建立约束方程，本算例是几何体积约束;
%%   4) 终止准则   结束迭代的条件，本算例是目标函数变化量<=0.010;
%%   5) 载荷工况   定义结构所有可能的受力情况，本算例是单一载荷;
%% ===================================================================================
%% 基于SIMP理论的优化准则法迭代分析流程:
%%   1) 定义设计域，选择合适的设计变量、目标函数以及约束函数等其他边界条件;
%%   2) 结构离散为有限元网格，计算优化前的单元刚度矩阵;
%%   3) 初始化单元设计变量，即给定设计域内的每个单元一个初始单元相对密度;
%%   4) 计算各离散单元的材料特性参数，计算单元刚度矩阵，组装结构总刚度矩阵，计算节点位移;
%%   5) 计算总体结构的柔度值及其敏度值，求解拉格朗日乘子;
%%   6) 用优化准则方法进行设计变量更新;
%%   7) 检查结果的收敛性，如未收敛则转 4)，循环迭代，如收敛则转 8);
%%   8) 输出目标函数值及设计变量值，结束计算。
%% ===================================================================================
 
function top99(nelx,nely,volfrac,penal,rmin);
% INITIALIZE
x(1:nely,1:nelx) = volfrac;      % 全局给一个初始体积分数，即初始相对密度
loop = 0;                        % 累计迭代次数，初始值为 0
change = 1.;                     % 该变量用来存储目标函数的变化量，如果小于某一个值就终止
% START ITERATION
while change > 0.001             % 改变量小于等于 0.01 时 (收敛)，跳出循环，否则执行
  loop = loop + 1;               % 每次循环会累计 1，也可以通过迭代多少次来决定程序的终止，而不仅仅是change
  xold = x;                      % xold 等于上一次循环的结果，用来和当前计算结果 x 比较，计算change
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);     % 根据单元数目和上次循环得到的 x，通过 f=ku 来计算 U (节点位移)，详见子函数
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;                     % ke 是单元的刚度矩阵，详见子函数（往下翻，在最后）
  c = 0.;                        % 用来叠加，先初始化为 0，这是目标函数
  % 遍历矩阵元素，先从左至右 (x 方向)，再从上至下 (y 方向), 
  % index 是 column by column from left to right, 
  % +  +  +
  % | /| /|
  % |/ |/ V
  % +  +  
  for ely = 1:nely               % y 方向遍历
    for elx = 1:nelx             % x 方向遍历
      n1 = (nely+1)*(elx-1)+ely; % (ely,elx) 单元的左上角节点编号
      n2 = (nely+1)* elx   +ely; % (ely,elx) 单元的右上角节点编号
                                 % 这里的编号，是从左到右，从上到下，比如最左上角为 1，往下 2，3，直到该列结束
      Ue = U([2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2],1); 
           % [LUx;    LUy ; RUx   ; RUy ; RDx   ; RDy   ; LDx   ; LDy   ];
           % 单元各节点编号的顺序, 各大写字母代表 L: Left, R: Right, U: Upper, D: Down, x y 分别代表 x y 方向的位移
           % 局部位移数组 Ue 储存 4 个节点共 8 个自由度位移，每个节点分别有 x、y 两个自由度; 
      c = c + x(ely,elx)^penal*Ue'*KE*Ue;
                                 % SIMP 模型，将设计变量 x 从离散型变成指函数型，指数就是惩罚因子 x(ely,elx)^penal
                                 % 通过每个单元的目标函数叠加得到总体的目标函数值，
                                 % 目标函数就是所有单元的柔度值，柔度就是 FU，其中 F=KU，所以总体来说就是 U'KU
                                 % 一般来说，大家都是以柔度作为目标函数，是因为 C=FU，
                                 % F 固定的情况下，结构肯定位移 U 越小越好；
                                 % U 固定的情况下，力越小越好，这样结构应力小，
                                 % 所以用 FU 的乘积能够比较好的衡量结构的刚度特性
                                 % c 是指 compliance
                                 % 这里目标函数指柔度 C=FU
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;  % 灵敏度，就是 c 对 x 求导，得到每 x 变化引起的 c 变化
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc); % 把刚计算的 dc 进行加权修正，防止棋盘格现象，具体参见后面子函数
                                       % 更新灵敏度
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD%利用优化准则更新变量
  [x]    = OC(nelx,nely,x,volfrac,dc); % 这一步是更新的关键部分，是这一部分让 x 朝着想要的方向 (优化准则) 更新，具体也参见后面子函数
                                       % 这里 x 指 xnew，省略 new 做连接，和 xold，实现迭代
% 不动点迭代法
% PRINT RESULTS
  change = max(max(abs(x-xold)));      % 正如上文提到的 xold，这里当前迭代中的 x 和之前一步的 xold 比较得到变化
% 如果结构基本没变化了，就可以认为迭代基本收敛了
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ... % 打印迭代的步数和目标函数
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...          % 打印体积分数，sum 加两次分别是行列相加
        ' ch.: ' sprintf('%6.3f',change )])                            % 打印每次更新后该变量的变化
% PLOT DENSITIES  
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6); % 这里画图，里面存在四种函数用法
  % colormap imagesc axis pause
  % colormap 查看并设置当前颜色图，gray 是灰度图，只有黑白两色，0-1 之间灰度变化，也有 winter（主要是蓝色，冷色系），summer 等色系，
  % 显示成灰度感觉最清晰
  % imagesc，imagesc(C) 将数组 C 中的数据显示为一个图像，该图像使用颜色图中的全部颜色。C 的每个元素指定图像的一个像素的颜色。
  % 生成的图像是一个 m×n 像素网格，其中 m 和 n 分别是 C 中的行数和列数。这些元素的行索引和列索引确定了对应像素的中心。
  % 举个例子，比如 0，1，2，4，7. 则颜色由黑到白（gray）分为 7 份，7 为黑，0 为白，1 2 4 对应相应的灰度程度按比例
  % axis 指坐标轴
  % pause(n) 暂停执行 n 秒，然后继续执行。必须启用暂停，此调用才能生效。
  % 只有暂停了，图像才能显示出来
end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 数学模型主要求解算法有：优化准则法(OC)、序列线性规划法(SLP)、移动渐进线法(MMA); 
% OC 适用于单约束优化问题求解，比如这里的"体积约束下的柔度最小化"问题，
% 当求解复杂的多约束拓扑优化问题时，采用 SLP 和 MMA 通常更方便
function [xnew]=OC(nelx,nely,x,volfrac,dc)  % xnew 是更新后的 x，xold 在程序初始已经储存了，需要参数 dc 和 volfrac
l1 = 0; l2 = 100000;                        % 利用二分法来逼近 l1 l2 边界
                                            % 这里我个人理解为类似于二分法，也不知道对不对，暂且这样理解吧
move = 0.2;                                 % 正向最大位移
while (l2-l1 > 1e-4)                        % 边界距离小于这个值时终止循环
  lmid = 0.5*(l2+l1);                       % 边界中间值
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid))))); 
  % -dc./lmid 为原文的 Be
  % 这个大小依次时这种情况，0.001 x-move xbnu x+move 1,如果满足这个大小顺序，就取 xbnu
  % 0.001 是最小密度下限，毕竟不能取负值，所以要大于这个值，故 max(0.001, )
  % x-move 是密度变化下限，每次不能变化太多，所以 max(m-move, )
  % 1 是上限，毕竟密度不能超过单位 1，所以 min( ,1)，这个满足了别急，还有一个上限
  % x+move 是变化的上限，所以 min( ,x+move)
  % 一切都是在判断 xbnu 的值，让它处在 x-move 和 x+move 之间，否则取上下限
  % 之所以这样做一方面是防止密度变成负值，或者变成大于 1 的值（0 无材料，1 实体材料，0-1 中间材料，超过该范围就没意义了）
  % 而 x-move and x+move should be included in [0.001,1]，因为有 x 比较小的情况，这时再 -0.2，比如 x-move = 0.13-0.2 = 负值
  % 负值就不行了，哪有密度是负的，所以取 0.001，1 也是同理
  if sum(sum(xnew)) - volfrac*nelx*nely > 0 % 最后判断体积，如果大于给定的体积分数，就取下限为中间值，整个区间往上逼近
                                            % 意思就是 limd 变大，而 limd 在分母上，整个体积就会往变小的方向移动
                                            % 否则小于给定的体积分数往下逼近
                                            % 意思就是如果体积小了，就增加体积，如果大，就减少
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)   % 网格独立性过滤器，其实就是检查灵敏度
                                            % 对过滤半径以内的其他单元的敏度 进行加权平均，修正单元的灵敏度，从而避免棋盘格现象，
                                            % 通常能够减少棋盘格现象的方法都能减少网格依赖性。
dcn=zeros(nely,nelx);                       % 先建立一个全部是 0 的用来待会存放的矩阵，用来存放灵敏度
                                            % 存放后方便使用 sum 符号，要进行累加来套公式
                                            % 提前规定好变量的空间，能提升效率
for i = 1:nelx
  for j = 1:nely                            % 对每一个单元进行循环
    sum=0.0;                                % 初始化总和，灵敏度中，每一个单元的Sum(H^)不同
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx) % 这两行表示，如果在圆的范围内，就要参与灵敏度修正
        % max(i-floor(rmin),1) 是为了防止超过边界，如果超过边界变成负值，就按 1 处理，min 同理，
        % floor 用来修正为整数
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely) % 注意 i,j 是目前所在的单元，k,l 是灵敏度所扫描的单元范围。
        fac = rmin-sqrt((i-k)^2+(j-l)^2); % fac 为 rmin - (ij到 kl 的距离)。卷积算子表示权重 H^
        sum = sum+max(0,fac);             % 用来表示权重则只有正值，若负值则看成 0，累加形成 sum 的效果 sum(H)，作分母的一部分
                                          % 理解上 max 无视，直接当成 fac 看就行
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k); % 这里累加形成分子，sum(H*x*dc/dx)
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);     % 这里将分子除以分母
    % Sigmund 提出的敏度过滤技术的本质是利用过滤范围内所有单元的敏度信息，修正中心单元的敏度信息，这里修正的思想不是太理解
    % 大概就是采用过滤半径范围内各单元敏度的加权平均值代替中心单元的敏度值
    % 为什么这种方法有效呢，你想想，其实棋盘格就是相邻两个单元一个想减少，一个想增大，我一加权，考虑局部整体
    % 变化就趋向于一致，就抑制了棋盘格
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)        % 这里是一个有限元程序，里面储存了边界信息
[KE] = lk;                                % 引入 ke
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)); % use "sparse" to save the memory maybe.
% generates an (m,n) all zero sparse matrix.
% here 2(nely+!)(+1) means
% 每个节点有两个方向的 k，就像是水平和竖直的弹簧，不知道这样理解对不对
% 如果以后用得到，再具体复习有限元这一部分内容
% 第二个 2(nely+!)(+1) 意思是其它节点对该节点处的影响，任何一个都可能在该处叠加，所以又一个 2(nely+!)(+1)，矩阵维数很大
F = sparse(2*(nely+1)*(nelx+1),1);  % 所有节点的力 (x,y 两个方向) 的列向量, 初始值为 0
U = zeros(2*(nely+1)*(nelx+1),1);   % 所有节点的位移 (x,y 两个方向) 列向量, 初始值为 0
% (nely+1)(nelx+1) 为总共节点数量
for elx = 1:nelx
  for ely = 1:nely                % 对单元循环
    n1 = (nely+1)*(elx-1)+ely;    % 定义一个单元上左坐标为n1，右上坐标为n2
    n2 = (nely+1)* elx   +ely;    % 先横向，再竖向，横向2n1-1，竖向2n1，按竖向增加，节点单元都是这样
                                  % 故左下为2n1+1 2n1+2
                                  % 用n1 n2可以简化下面代码，不然看起来太长
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2]; % 这里以左上为起点，顺时针排列，感觉应该不影响
                                                                         % 顺时针逆时针都一样，数值不变，都是对应的
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE; % 在对应位置处累加，注意这里 edof 中位置不一定相邻，
                                                       % 这里是在对单元的刚度矩阵，对号放入到总体的刚度矩阵中
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F(2,1) = -1;                      % 负载设定，-1 表示 y 负方向单位力，具体图见 Sigmund 文献资料
                                  % (index, 列)
fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); % 这里的 union 按从小到大返回并集，
                                                             % 1:2:2*(nely+1) 表示所有第一列的 x，
                                                             % 2*(nelx+1)*(nely+1) 表示右下角节点 y
                                                             % 表示左边界 x 固定，右下角 y 固定
alldofs     = [1:2*(nely+1)*(nelx+1)];                       % 列举所有的自由度 (index)
freedofs    = setdiff(alldofs,fixeddofs);                    % 除去固定自由度，剩下的完全自由的节点自由度 (index)
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);        % 求解自由的自由度，有限元中已经练习过，可以整体拿掉某行某列进行计算
U(fixeddofs,:)= 0;                                           % 固定处的位移为0
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk                                             % 这里就是储存了一个 ke 的信息，这里采用了平面 4 节点的有限元单元
                                                             % ke 为单元的刚度矩阵，ke x u = p, (u 为节点位移，p 为节点外力)
E = 1.;                                                      % 材料杨氏弹性模量
nu = 0.3;                                                    % 材料 poisson's ratio
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...           % 4 节点共有 8 个 DOF
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
% 有限元方法计算的一个重要的系数矩阵，表征单元体的受力与变形关系;
% 特点：对称性、奇异性、主对角元素恒正、奇数（偶数）行和为0;
% 矩形单元 4 节点 8*8矩阵;
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)    % 直接套公式就行
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
