close
clc
clear
xdata=linspace(0,2*pi,15);
y=5*sin(xdata)+2*xdata+xdata.^2;
y=y+2*rand(1,15);
plot(xdata,y,'o') %���Ƴ�ԭʼ���ݵ�ɢ��ͼ
hold on
 
fun=@(x,xdata) x(1)*sin(xdata)+x(2)*xdata+x(3)*xdata.^2;%������Ϻ�������ʽ
x=lsqcurvefit(fun,[0 0 0],xdata,y);% [0 0 0]Ϊ��ֵ����ʽ��ʼϵ��a0��
%xdataΪ����ԭʼ����x���꣬yΪԭʼ���������꣬����ֵxΪ����ϵ������
xx=linspace(0,2*pi,150);%��Ϻ������ߵ�x����
yy=fun(x,xx);%��Ϻ�������
plot(xx,yy)%������Ϻ�������
 
lb=[-1 -1 -1];%ϵ������
ub=[6 3 2]; %ϵ������
x=lsqcurvefit(fun,[0 0 0],xdata,y,lb,ub);
xx=linspace(0,2*pi,150);
yy=fun(x,xx);
plot(xx,yy)
%ѡ������
%options = optimset(Name,Value,...)�����Խ������������Ϊ��һ������
%options = optimoptions(SolverName,Name,Value,...)���Խ������������Ϊ��һ������
% options=optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');%�㷨����
options=optimoptions('lsqcurvefit','Display','final');%��ʾ����
lb=[-1 -1 -1];
ub=[6 3 2];
[x,~,~,exitflag,~,~,jacobian]=lsqcurvefit(fun,[0 0 0],xdata,y,lb,ub,options);
xx=linspace(0,2*pi,150);
yy=fun(x,xx);
plot(xx,yy)