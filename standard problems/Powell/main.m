%线圈的信息
x=0;            % 当前线圈的x坐标
y=0;            % 当前线圈的y坐标
h=0;
CoilPostion=[x y h]; % 线圈位置，默认高度为0   
% theta_c = unidrnd(30)*pi/180;
% phi_c = unidrnd(30)*pi/180;
I=20;
R=0.4;
f=1000;

%目标物体的信息
CR=0.1;
p=0.8/(2*CR);
Cu=1;
Csigma=5.71*10^7;
postion=[0 0 -2];%目标位置 [-2 2 -1]
thetatrue=20;
phitrue=30;
theta=thetatrue*pi/180;
phi=phitrue*pi/180;

v_x=-2:0.5:2;
v_y=-2:0.5:2;

[m_x, m_y]=meshgrid(v_x, v_y);
Hx=zeros(size(m_x));
Hy=zeros(size(m_x));
Hz=zeros(size(m_x));
for i=1:length(v_y)          % 对一个探测点的线圈姿态进行变化，得到一组探测数据  
    for j=1:length(v_x)
        x=v_x(j);            % 当前线圈的x坐标
        y=v_y(i);            % 当前线圈的y坐标
        CoilPostion=[x y h];
%         theta_c = unidrnd(30)*pi/180;
%         phi_c = unidrnd(30)*pi/180;
        
        [CMxyz,M]=MomCylinder(FirstField(I,R,f,(postion-CoilPostion)),Cu,Csigma,f,p,theta,phi,CR); % 重复出现多行数据的问题语句 
        CMM=[postion(:) CMxyz(:)];
        [Hx(i,j),Hy(i,j),Hz(i,j)]=HFieldModel(CMM,x,y,h);
    end
end
v_M=[M(1,1);M(2,2);M(3,3);M(1,2);M(1,3);M(2,3);];  % 取张量矩阵中的6个独立元素
AlphaTrue=[postion(:)' M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3) ]; % 真实值

[betaTrue,~]=GetBetaAndAngle(v_M);  % 获取探测目标的主轴极化率
            
% nHx=awgn(abs(Hx),30,'measured');
% nHy=awgn(abs(Hy),30,'measured');
% nHz=awgn(abs(Hz),30,'measured');
nHx=awgn((Hx),10,'measured');
nHy=awgn((Hy),10,'measured');
nHz=awgn((Hz),10,'measured');
k=0;
m_pk=zeros(length(v_x)*length(v_y),3);   % 线圈位置集合
m_data=zeros(length(v_x)*length(v_y),3); % 观测数据集
% m_coils=zeros(length(v_x)*length(v_y),2);
for i=1:length(v_y)          % 对二维区域的单元格进行扫描，相对于移动线圈
    for j=1:length(v_x)
        k=k+1;
        x=v_x(j);            % 当前线圈的x坐标
        y=v_y(i);            % 当前线圈的y坐标
        m_data(k,:)=[nHx(i,j) nHy(i,j) nHz(i,j)];    % 将观测数据从矩阵转换成一数组
        m_pk(k,:)=[x y 0];                                              % 观测数据对应的探测器位置
%         m_coils(k,:)=[coils_theta(i,j) coils_phi(i,j)];
    end
end
m_pk=m_pk(1:k,:);
m_data=m_data(1:k,:);
% m_coils = m_coils(1:k,:);
                
v_initial=[0 0 -5 0 0 0 0 0 0]; %1.换初值 2.看一下出现奇异值点的点 有没有共同点 
% fun = @(x)ObjectFun3_1(x,m_pk,m_data);
v_initial=[-1.2 1];
fun = @(x) [100*(x(2) - x(1)^2)^2,(1 - x(1))^2]; % Rosenbrock
fprime = @(x)GetObjGrad3_1(x,m_pk,m_data);
t1=clock;
%[n,x] = ConjugateGradient(fun,fprime,v_initial,1000,1e-9);
%[n,x] = GradientDescent(fun,fprime,1000,v_initial,5000,1e-9);
% [n,x] = QuasiNewtonBFGS(fun,fprime,v_initial,1000,1e-9);
%[n,x] = Newton(fun,fprime,hess,v_initial,1000,1e-9);
%[n,x] = SteepestDescent(fun,fprime,v_initial,10000,1e-9);   
% [n,x] = LevenbergMarquardt_1(fun,fprime,v_initial,1000,1e-9,m_pk); % 异常输出语句
[n,x] = LevenbergMarquardt_0(fun,fprime,v_initial,1000,1e-9,m_pk); % 异常输出语句 返回迭代次数和值
t2=clock;
runtime = etime(t2,t1);
%% 结果处理
% Scene.result.Alpha=x;
% Scene.result.Alpha0=v_initial;
e_cf=norm(GetResidualVector_1(x,m_pk,m_data),2);
i_cf=norm(GetResidualVector_1(v_initial,m_pk,m_data),2);
t_cf=norm(GetResidualVector_1(AlphaTrue,m_pk,m_data),2);
e_pos=x(1:3); % 定位结果
[betas,angles]=GetBetaAndAngle(x(4:9));  % 姿态和极化率估计结果
pos_e=[abs(e_pos(1)-AlphaTrue(1)) abs(e_pos(2)-AlphaTrue(2)) abs(e_pos(3)-AlphaTrue(3))];
betax = abs(betas(1)-betaTrue(1));
% betaTrue(1)
% betaTrue(2)
% betaTrue(3)
% betas(1)
% betas(2)
% betas(3)
betay = abs(betas(2)-betaTrue(2));
betaz = abs(betas(3)-betaTrue(3));
theta=abs(angles(1)-thetatrue);
phi=abs(angles(2)-phitrue);
psi=0;

shuchu = ('【lm结果】');
disp(shuchu);
shuchu=['迭代次数：',num2str(n),' 运行时间：',num2str(runtime)];
disp(shuchu);
shuchu=['迭代初值：', num2str(v_initial)];
disp(shuchu);
shuchu=['迭代结果：', num2str(x)];
disp(shuchu);
shuchu=['真实值：', num2str(AlphaTrue)];
disp(shuchu);
shuchu=['误差：', num2str(abs(AlphaTrue-x))];
disp(shuchu);
shuchu=['【目标函数值】'];
disp(shuchu);
shuchu=['迭代初值处：', num2str(i_cf)];
disp(shuchu);
shuchu=['收敛值处：', num2str(e_cf)];
disp(shuchu);
shuchu=['真实值处：', num2str(t_cf)];
disp(shuchu);
shuchu=['【定位误差/m】'];
disp(shuchu);
shuchu=['x,y,z=',  num2str(pos_e)];
disp(shuchu);
shuchu=['【主轴极化率误差/百分比】'];
disp(shuchu);
shuchu=['x,y,z=',  num2str(betax), ' ' ,num2str(betay), ' ', num2str(betaz)];
disp(shuchu);
shuchu=['【姿态误差/角度】'];
disp(shuchu);
shuchu=['theta,phi,psi=',  num2str(theta), ' ' ,num2str(phi), ' ', num2str(psi)];
disp(shuchu);

% filename1 = 'C:\Users\thinkpad\Desktop\1\毕设\实验.csv';
% fid1 = fopen(filename1, 'at+');
% fprintf(fid1,['%d,'],n); %频率
% fprintf(fid1,['%d,'],runtime);
% fprintf(fid1,['%d,'],pos_e(1));
% fprintf(fid1,['%d,'],pos_e(2));
% fprintf(fid1,['%d,'],pos_e(3));
% fprintf(fid1,['%d,'],betax);
% fprintf(fid1,['%d,'],betay);
% fprintf(fid1,['%d,'],betaz);
% fprintf(fid1,['%d,'],theta);
% fprintf(fid1,['%d,'],phi);
% fprintf(fid1,['%d,'],i_cf);
% fprintf(fid1,['%d,'],e_cf);
% fprintf(fid1,['%d,'],t_cf);
% fprintf(fid1,['%d,'],AlphaTrue(1));
% fprintf(fid1,['%d,'],AlphaTrue(2));
% fprintf(fid1,['%d,'],AlphaTrue(3));
% fprintf(fid1,['%d,'],betaTrue(1));
% fprintf(fid1,['%d,'],betaTrue(2));
% fprintf(fid1,['%d,'],betaTrue(3));
% fprintf(fid1,['%d,'],thetatrue);
% fprintf(fid1,['%d\n'],phitrue);
% std=fclose('all');     
