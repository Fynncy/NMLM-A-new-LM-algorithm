%��Ȧ����Ϣ
x=0;            % ��ǰ��Ȧ��x����
y=0;            % ��ǰ��Ȧ��y����
h=0;
CoilPostion=[x y h]; % ��Ȧλ�ã�Ĭ�ϸ߶�Ϊ0   
% theta_c = unidrnd(30)*pi/180;
% phi_c = unidrnd(30)*pi/180;
I=20;
R=0.4;
f=1000;

%Ŀ���������Ϣ
CR=0.1;
p=0.8/(2*CR);
Cu=1;
Csigma=5.71*10^7;
postion=[0 0 -2];%Ŀ��λ�� [-2 2 -1]
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
for i=1:length(v_y)          % ��һ��̽������Ȧ��̬���б仯���õ�һ��̽������  
    for j=1:length(v_x)
        x=v_x(j);            % ��ǰ��Ȧ��x����
        y=v_y(i);            % ��ǰ��Ȧ��y����
        CoilPostion=[x y h];
%         theta_c = unidrnd(30)*pi/180;
%         phi_c = unidrnd(30)*pi/180;
        
        [CMxyz,M]=MomCylinder(FirstField(I,R,f,(postion-CoilPostion)),Cu,Csigma,f,p,theta,phi,CR); % �ظ����ֶ������ݵ�������� 
        CMM=[postion(:) CMxyz(:)];
        [Hx(i,j),Hy(i,j),Hz(i,j)]=HFieldModel(CMM,x,y,h);
    end
end
v_M=[M(1,1);M(2,2);M(3,3);M(1,2);M(1,3);M(2,3);];  % ȡ���������е�6������Ԫ��
AlphaTrue=[postion(:)' M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3) ]; % ��ʵֵ

[betaTrue,~]=GetBetaAndAngle(v_M);  % ��ȡ̽��Ŀ������Ἣ����
            
% nHx=awgn(abs(Hx),30,'measured');
% nHy=awgn(abs(Hy),30,'measured');
% nHz=awgn(abs(Hz),30,'measured');
nHx=awgn((Hx),10,'measured');
nHy=awgn((Hy),10,'measured');
nHz=awgn((Hz),10,'measured');
k=0;
m_pk=zeros(length(v_x)*length(v_y),3);   % ��Ȧλ�ü���
m_data=zeros(length(v_x)*length(v_y),3); % �۲����ݼ�
% m_coils=zeros(length(v_x)*length(v_y),2);
for i=1:length(v_y)          % �Զ�ά����ĵ�Ԫ�����ɨ�裬������ƶ���Ȧ
    for j=1:length(v_x)
        k=k+1;
        x=v_x(j);            % ��ǰ��Ȧ��x����
        y=v_y(i);            % ��ǰ��Ȧ��y����
        m_data(k,:)=[nHx(i,j) nHy(i,j) nHz(i,j)];    % ���۲����ݴӾ���ת����һ����
        m_pk(k,:)=[x y 0];                                              % �۲����ݶ�Ӧ��̽����λ��
%         m_coils(k,:)=[coils_theta(i,j) coils_phi(i,j)];
    end
end
m_pk=m_pk(1:k,:);
m_data=m_data(1:k,:);
% m_coils = m_coils(1:k,:);
                
v_initial=[0 0 -5 0 0 0 0 0 0]; %1.����ֵ 2.��һ�³�������ֵ��ĵ� ��û�й�ͬ�� 
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
% [n,x] = LevenbergMarquardt_1(fun,fprime,v_initial,1000,1e-9,m_pk); % �쳣������
[n,x] = LevenbergMarquardt_0(fun,fprime,v_initial,1000,1e-9,m_pk); % �쳣������ ���ص���������ֵ
t2=clock;
runtime = etime(t2,t1);
%% �������
% Scene.result.Alpha=x;
% Scene.result.Alpha0=v_initial;
e_cf=norm(GetResidualVector_1(x,m_pk,m_data),2);
i_cf=norm(GetResidualVector_1(v_initial,m_pk,m_data),2);
t_cf=norm(GetResidualVector_1(AlphaTrue,m_pk,m_data),2);
e_pos=x(1:3); % ��λ���
[betas,angles]=GetBetaAndAngle(x(4:9));  % ��̬�ͼ����ʹ��ƽ��
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

shuchu = ('��lm�����');
disp(shuchu);
shuchu=['����������',num2str(n),' ����ʱ�䣺',num2str(runtime)];
disp(shuchu);
shuchu=['������ֵ��', num2str(v_initial)];
disp(shuchu);
shuchu=['���������', num2str(x)];
disp(shuchu);
shuchu=['��ʵֵ��', num2str(AlphaTrue)];
disp(shuchu);
shuchu=['��', num2str(abs(AlphaTrue-x))];
disp(shuchu);
shuchu=['��Ŀ�꺯��ֵ��'];
disp(shuchu);
shuchu=['������ֵ����', num2str(i_cf)];
disp(shuchu);
shuchu=['����ֵ����', num2str(e_cf)];
disp(shuchu);
shuchu=['��ʵֵ����', num2str(t_cf)];
disp(shuchu);
shuchu=['����λ���/m��'];
disp(shuchu);
shuchu=['x,y,z=',  num2str(pos_e)];
disp(shuchu);
shuchu=['�����Ἣ�������/�ٷֱȡ�'];
disp(shuchu);
shuchu=['x,y,z=',  num2str(betax), ' ' ,num2str(betay), ' ', num2str(betaz)];
disp(shuchu);
shuchu=['����̬���/�Ƕȡ�'];
disp(shuchu);
shuchu=['theta,phi,psi=',  num2str(theta), ' ' ,num2str(phi), ' ', num2str(psi)];
disp(shuchu);

% filename1 = 'C:\Users\thinkpad\Desktop\1\����\ʵ��.csv';
% fid1 = fopen(filename1, 'at+');
% fprintf(fid1,['%d,'],n); %Ƶ��
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
