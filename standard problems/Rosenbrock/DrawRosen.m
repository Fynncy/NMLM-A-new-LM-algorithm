function DrawRosen(xlog, ylog, xdlog, n)
% %UNTITLED5 此处显示有关此函数的摘要
% %   此处显示详细说明
% 
s = 0.02;
x1 = [-2 : s : 2+s];
x2 = [-1 : s : 3+s];
[x1, x2] = meshgrid(x1, x2);
y = 100 .* (x2 - x1.^2).^2 + (1 - x1).^2;
miny = min(y(:));
maxy = max(y(:));
C = miny + (maxy-miny).*log(1+y-miny)./log(1+maxy-miny);

figure;
surf(x1, x2, y, C, 'EdgeColor', 'none', 'LineStyle', 'none');
axis([-2, 2, -1, 3, 0, 2500]);%坐标轴的范围
xlabel('x1');ylabel('x2');zlabel('y','Rotation',0);
title([name ' optimizing track 3D']);
%colormap jet;

hold on
plot3(xlog(1,:),xlog(2,:),ylog,'-ro','LineWidth',1.5,'MarkerSize',3,'MarkerFaceColor','r');
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
plot(xlog(1,:),xlog(2,:),'-r*','LineWidth',1.5,'MarkerSize',3,'MarkerFaceColor','r');
title([name ' optimizing track in contour']);
xlabel('x1');ylabel('x2','Rotation',0);
grid on
% grid minor
% 
figure;
plot(1:n, ylog,'k','LineWidth',2);
axis([-inf,inf, -inf,inf]);
title([name ' y value variations with the iterations']);
xlabel('Number of iterations');ylabel('y','Rotation',0);
grid on
% set(gca, 'GridLineStyle', '--'); 
% 
figure;
plot(1:n, xdlog,'k','LineWidth',2);
axis([-inf,inf, -inf,inf]);
title([name ' x diff value variations with the iterations']);
xlabel('Number of iterations');ylabel('x diff');
grid on

% figure;
% plot(xlog(1,:),ylog,'g',xlog(2,:),ylog,'b--o');


% % 绘制随时间变化情况
% t = seconds(0:0.5);
% y = 1:5;
% plot(t,y,'DurationTickFormat','hh:mm:ss')


end

