function DrawPowell(xlog, ylog, xdlog, n, name)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明

s = 0.02;
x1 = [-2 : s : 3+s];
x2 = [-2 : s : 3+s];
x3 = [-1 : s : 4+s];
x4 = [-3 : s : 2+s];
[x1, x2] = meshgrid(x1, x2);
% y = 100 .* (x2 - x1.^2).^2 + (1 - x1).^2;
y = (x1 + 10.* x2).^2 + 5 .* (x3 - x4).^2 + (x2 - 2 .* x3).^4 + 10 .* (x1 - x4).^4;
miny = min(y(:));
maxy = max(y(:));
C = miny + (maxy-miny).*log(1+y-miny)./log(1+maxy-miny); %用于指定在不同高度下的颜色范围

% figure(101);  % plot convergence history of parameters, reduced chi^2, lambda
%  clf
% subplot(221) %subplot(x,y,n)
%     surf(x1, x2, y, C, 'EdgeColor', 'none', 'LineStyle', 'none');
% subplot(222)
%     surf(x1, x3, y, C, 'EdgeColor', 'none', 'LineStyle', 'none');
% subplot(223)
%     surf(x2, x3, y, C, 'EdgeColor', 'none', 'LineStyle', 'none');
% subplot(224)
%     surf(x3, x4, y, C, 'EdgeColor', 'none', 'LineStyle', 'none');




figure;
% surf(x1, x2, y, C, 'EdgeColor', 'none', 'LineStyle', 'none'); %x,y,z是同维矩阵，x,y是网格坐标矩阵，z是网格点的高度矩阵
% % 当x,y是向量时，要求x的长度必须等于z矩阵的列数，y的长度等于z矩阵的行数，
% % x，y向量元素的组合构成网格点的x，y坐标，z坐标则取自z矩阵，然后绘制三维曲面图。
% axis([-1, 4, -1, 4, 0, 2000]);%坐标轴的范围
% xlabel('x1');ylabel('x2');zlabel('y','Rotation',0);
% title([name ' optimizing track 3D']);
% %colormap jet;

hold on
plot3(xlog(1,:),xlog(2,:),ylog,'-ro','LineWidth',2,'MarkerSize',3,'MarkerFaceColor','r');
colorbar


%legend(temp,{'track'},'Location','northwest');
% 
figure;
%[~, h] = contourf(x1, x2, C, 50);
%set(h,'LineStyle','none')
%alpha(gca,0.8);
%hold on;
[~, h] = contour(x1, x2, y);
set(h,'ShowText','on','LineWidth',1,'LineColor','k','LevelList',[ 0.1 1 2.5 10 100 1000],'TextList',[1 10 100 1000])

hold on
plot(xlog(1,:),xlog(4,:),'-r*','LineWidth',1.5,'MarkerSize',3,'MarkerFaceColor','r');
title([name ' optimizing track in contour']);
xlabel('x1');ylabel('x2','Rotation',0);
grid on
% grid minor
% 
% figure;
% plot(1:n, ylog,'k','LineWidth',2);
% axis([-inf,inf, -inf,inf]);
% title([name ' y value variations with the iterations']);
% xlabel('Number of iterations');ylabel('y','Rotation',0);
% grid on
% % set(gca, 'GridLineStyle', '--'); 
% % 
% figure;
% plot(1:n, xdlog,'k','LineWidth',2);
% axis([-inf,inf, -inf,inf]);
% title([name ' x diff value variations with the iterations']);
% xlabel('Number of iterations');ylabel('x diff');
% grid on

figure;
% plot(xlog(1,:),ylog,'g',xlog(2,:),ylog,'b--o');
plot(1:n, xlog(1,:),'g--o', 1:n, xlog(2,:),'b',1:n, xlog(3,:),'r',1:n, xlog(4,:),'y');


figure(101);  % plot convergence history of parameters, reduced chi^2, lambda
 clf
subplot(221) %subplot(x,y,n)
    plot(1:n, xlog(1,:),'g','LineWidth',2);
subplot(222)
    plot(1:n, xlog(2,:),'b','LineWidth',2);
subplot(223)
    plot(1:n, xlog(3,:),'r','LineWidth',2);
subplot(224)
    plot(1:n, xlog(4,:),'y','LineWidth',2);


end

