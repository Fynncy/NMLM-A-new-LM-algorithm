function [grad] = JacobianPowell(x)
%JacobianRosen 给定包含两个元素的一位向量x，返回Rosenbrock在x点的梯度grad
% 分别对x1和x2求导
%   此处显示详细说明

% grad = [200*(x(2)-x(1).^2)*(-2*x(1)) + 2*(x(1)-1), 200*(x(2)-x(1).^2)];
grad = [2 * (x(1)+ 10 * x(2)) + 40 * (x(1) - x(4)).^3, 20 * (x(1) + 10 * x(2)) + 4 * (x(2) - 2 * x(3)).^3, ...
        10 * (x(3) - x(4)) - 8 * (x(2) - 2 * x(3)).^3, -10 * (x(3) - x(4)) - 40 * (x(1) - x(4)).^3];
end

