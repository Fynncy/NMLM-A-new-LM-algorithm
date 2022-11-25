clc;
clear;
close all;

%建议阅读顺序：main-Rosenbrock-JacobianRosen-HessianRosen-GradientDescent-
%SteepestDescent-StepLength-Newton-QuasiNewtonGFGS-ConjugateGradient.

%定义函数，这相当于将函数指针赋值给f，Rosenbrock是在别的文件中写好的函数
f = @Powell;
%梯度
grad = @JacobianPowell;
%海森矩阵
hessian = @HessianPowell;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rosenbrock Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=round(rand(1,1)*10)-5;
b=round(rand(1,1)*5)-2;
c=round(rand(1,1)*4)-1;
d=round(rand(1,1)*6)-2;

% x1=-5:1:4;
% x2=-2:1:3;
% x1=-2:0.1:2;
% x2=-2:0.1:2;
x0 = [a, b, c, d];

% [xlog, ylog, xdlog, n] = LevenbergMarquardt_m1(f, grad, [3,-1,0,1], 100, 1e-4, hessian);
t1=clock;
[xlog, ylog, xdlog, n] = LevenbergMarquardt_m1(f, grad, [0,0,0,0], 100, 1e-4, hessian);
t2=clock;
fprintf("time=%f\n",etime(t2,t1));
% fprintf(fid1,['%s,'],method);
x = xlog';
disp('======================================================');
disp(x0);
disp('======================================================');
disp(x(n,:));
disp('======================================================');
disp(ylog(1,n));
disp('======================================================');
% disp(xdlog);
% disp('======================================================');
disp(n);
% B=reshape(xlog,[],2);
% disp(B);
% DrawPowell(xlog, ylog, xdlog, n, 'Powell singular function');
lylog = length(ylog);
row = zeros(1,length(ylog));
% row = [0,0,0,0,0,0,0,0,0,0,0,0,0];
ylog_com = [ylog;row];
ylog_com = [ylog_com;row];
ylog_com = [ylog_com;row];
final = [xlog',ylog_com'];
filename1 = 'D:\博士项目研究\目标探测\论文\优化\反演\LM算法\matlab_result\final_result\Powell\result.csv';
fopen(filename1, 'at+');
fid1=fopen(filename1, 'at+');
% fprintf(fid1,['%s,'],'xlog');
dlmwrite(filename1, x0);
fprintf(fid1,['%s\n,'],'');
dlmwrite(filename1, final);
% dlmwrite(filename1, xlog');
% fprintf(fid1,['%s\n,'],'');
% fprintf(fid1,['%s,'],'ylog');
% dlmwrite(filename1, ylog);
std=fclose('all'); 

fprintf('finish');

% x1=-2:0.1:2;
% x2=-2:0.1:2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rosenbrock Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% [xlog, ylog, xdlog, n] = LevenbergMarquardt_m1(f, grad, [-1.2,1], 100, 1e-4, hessian);
% toc
% DrawRosen(xlog, ylog, xdlog, n, 'Rosenbrock banana function');

%%%%%时间和迭代次数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % tic;
% filename1 = 'D:\博士项目研究\目标探测\论文\优化\反演\LM算法\matlab_result\final_result\Rosenbrock\result.csv';
% fid1 = fopen(filename1, 'at+');
% % fprintf(fid1,['%s,'],'mLM1');
% for i=1:50
%     x1=-100:0.1:5;
%     x2=-2:0.1:20;
%     a=randi(length(x1));
%     b=randi(length(x2));
%     t1 = clock;
%     
% %     [xlog, ylog, xdlog, n] = LevenbergMarquardt_m1(f, grad, [x1(a),x2(b)], 1000, 1e-4, hessian);
%     [xlog, ylog, xdlog, n] = LevenbergMarquardt_m1(f, grad, [5,20], 1000, 1e-4, hessian);
%     t2 = clock;
%     fprintf("num=%d, x1=%f, x2=%f, time=%f, iter=%d\n",i, x1(a), x2(b), etime(t2,t1), n);
% %     fprintf(fid1,['%f,'],x1(a));
% %     fprintf(fid1,['%f,'],x2(b));
% %     fprintf(fid1,['%f,'],etime(t2,t1));
% %     fprintf(fid1,['%d\n'],n);
% end
% std=fclose('all');  
% % toc
% % DrawRosen(xlog, ylog, xdlog, n, 'Rosenbrock banana function');
%%%%%%%%时间和迭代次数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vfun = @(x) [100*(x(2) - x(1)^2)^2 + (1 - x(1))^2];
% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
% tic
% [x] = lsqnonlin(vfun,[2,3],[],[],options);
% toc

% x
% unloadlibrary fortran

% figure;
% f = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
% x = linspace(-1.5,1.5); y = linspace(-1,3);
% [xx,yy] = meshgrid(x,y); ff = f(xx,yy);
% levels = 10:10:300;
% LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
% figure, contour(x,y,ff,levels,LW,1.2), colorbar
% axis([-1.5 1.5 -1 3]), axis square, hold on

% % 绘制拟合曲线
% randn('seed',0);
% consts = [ ];
% Npnt = 100;				  % number of data points
% t = [1:Npnt]';				  % independent variable, column vector
% y_dat = f([2,3]);
% y_dat = y_dat + 0.5*randn(Npnt,1);
% [xlog, ylog, xdlog, n, xk] = LevenbergMarquardt_m1(f, grad, [-1.2,1], 100, 1e-4, hessian);
% y_fit = f(xk);
% sigma_y=[];
% for i=1:Npnt
%         sigma_y(i)=i*0.01*randi(10);
% end
% 
% figure(101); % ------------ plot data, fit, and confidence interval of fit
% confidence_level = 0.99;   % confidence level for error confidence interval;
% z = norminv((1+confidence_level)/2);
% clf
%    plot(t,y_dat,'og', t,y_fit,'-b', t,y_fit+z*sigma_y,'.k', t,y_fit-z*sigma_y,'.k');
% % % subplot(411)
% % %    plot(t,y_dat,'og');
% % %    ylabel('y_dat')
% % % subplot(412)
% % %    plot(t,y_fit,'-b');
% % %    ylabel('y_fit')
% % % subplot(413)
% % %    plot(t,y_fit+z*sigma_y,'.k');
% % %    ylabel('sandian')
% % % subplot(414)
% % %    plot(t,y_fit-z*sigma_y,'.k');
% % %    ylabel('sandian')
% % legend('y_{data}','y_{fit}','99% c.i.','',0);
% % ylabel('y(t)')
% xlabel('t')





