close
clc
clear
xdata=linspace(0,2*pi,15);
y=5*sin(xdata)+2*xdata+xdata.^2;
y=y+2*rand(1,15);
plot(xdata,y,'o') %绘制出原始数据的散点图
hold on
 
fun=@(x,xdata) x(1)*sin(xdata)+x(2)*xdata+x(3)*xdata.^2;%待求拟合函数的形式
x=lsqcurvefit(fun,[0 0 0],xdata,y);% [0 0 0]为插值多项式初始系数a0，
%xdata为输入原始数据x坐标，y为原始数据纵坐标，返回值x为待求系数矩阵
xx=linspace(0,2*pi,150);%拟合函数曲线的x坐标
yy=fun(x,xx);%拟合函数曲线
plot(xx,yy)%绘制拟合函数曲线
 
lb=[-1 -1 -1];%系数下限
ub=[6 3 2]; %系数上限
x=lsqcurvefit(fun,[0 0 0],xdata,y,lb,ub);
xx=linspace(0,2*pi,150);
yy=fun(x,xx);
plot(xx,yy)
%选项设置
%options = optimset(Name,Value,...)不可以将求解器名称作为第一个参数
%options = optimoptions(SolverName,Name,Value,...)可以将求解器名称作为第一个参数
% options=optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');%算法设置
options=optimoptions('lsqcurvefit','Display','final');%显示设置
lb=[-1 -1 -1];
ub=[6 3 2];
[x,~,~,exitflag,~,~,jacobian]=lsqcurvefit(fun,[0 0 0],xdata,y,lb,ub,options);
xx=linspace(0,2*pi,150);
yy=fun(x,xx);
plot(xx,yy)