% % Generate some data.
% 
% data = randn(10000, 2);
% 
% % Scale and rotate the data (for demonstration purposes).
% 
% data(:,1) = data(:,1) * 2;
% 
% theta = deg2rad(130);
% 
% data = ([cos(theta) -sin(theta); sin(theta) cos(theta)] * data')';
% 
% % Get some info.
% 
% m = mean(data);
% 
% s = std(data);
% 
% axisMin = m - 4 * s;
% 
% axisMax = m + 4 * s;
% 
% % Plot data points on (X=data(x), Y=data(y), Z=0)
% 
% plot3(data(:,1), data(:,2), zeros(size(data,1),1), 'k.', 'MarkerSize', 1);
% 
% % Turn on hold to allow subsequent plots.
% 
% hold on
% 
% % Plot the ellipse using Eigenvectors and Eigenvalues.
% 
% data_zeroMean = bsxfun(@minus, data, m);
% 
% [V,D] = eig(data_zeroMean' * data_zeroMean / (size(data_zeroMean, 1)));
% 
% [D, order] = sort(diag(D), 'descend');
% 
% D = diag(D);
% 
% V = V(:, order);
% 
% V = V * sqrt(D);
% 
% t = linspace(0, 2 * pi);
% 
% e = bsxfun(@plus, 2*V * [cos(t); sin(t)], m');
% 
% plot3(e(1,:), e(2,:), zeros(1, 100), 'g-', 'LineWidth', 2); % xy平面上的圈
% 
% maxP = 0;
% 
% for side = 1:2
% 
% % Calculate the histogram.
% 
% p = [0 hist(data(:,side), 20) 0];
% 
% p = p / sum(p);
% 
% maxP = max([maxP p]);
% 
% dx = (axisMax(side) - axisMin(side)) / numel(p) / 2.3;
% 
% p2 = [zeros(1,numel(p)); p; p; zeros(1,numel(p))]; p2 = p2(:);
% 
% x = linspace(axisMin(side), axisMax(side), numel(p));
% 
% x2 = [x-dx; x-dx; x+dx; x+dx]; x2 = max(min(x2(:), axisMax(side)), axisMin(side));
% 
% % Calculate the curve.
% 
% nPtsCurve = numel(p) * 10;
% 
% xx = linspace(axisMin(side), axisMax(side), nPtsCurve);
% 
% % Plot the curve and the histogram.
% 
% if side == 1
% 
% plot3(xx, ones(1, nPtsCurve) * axisMax(3 - side), spline(x,p,xx), 'r-', 'LineWidth', 2);
% 
% plot3(x2, ones(numel(p2), 1) * axisMax(3 - side), p2, 'k-', 'LineWidth', 1);
% 
% else
% 
% plot3(ones(1, nPtsCurve) * axisMax(3 - side), xx, spline(x,p,xx), 'b-', 'LineWidth', 2);
% 
% plot3(ones(numel(p2), 1) * axisMax(3 - side), x2, p2, 'k-', 'LineWidth', 1);
% 
% end
% 
% end
% 
% % Turn off hold.
% 
% hold off
% 
% % Axis labels.
% 
% xlabel('x');
% 
% ylabel('y');
% 
% zlabel('p(.)');
% 
% axis([axisMin(1) axisMax(1) axisMin(2) axisMax(2) 0 maxP * 1.05]);
% 
% grid on;

%=========================================================================
% figure;
% [X,Y,Z] = peaks(25);
% CO(:,:,1) = zeros(25); % red
% CO(:,:,2) = ones(25).*linspace(0.5,0.6,25); % green
% CO(:,:,3) = ones(25).*linspace(0,1,25); % blue
% surf(X,Y,Z,CO);
% 
% figure;
% [X,Y] = meshgrid(-5:.5:5);
% Z = Y.*sin(X) - X.*cos(Y);
% s = surf(X,Y,Z,'FaceAlpha',0.5);
%=========================================================================


%降维
A=[100 200 300 400 500];
[A1,PS]=mapminmax(A);
disp(A);