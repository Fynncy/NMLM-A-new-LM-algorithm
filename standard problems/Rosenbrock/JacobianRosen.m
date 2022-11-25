function [grad] = JacobianRosen(x)
%JacobianRosen 给定包含两个元素的一位向量x，返回Rosenbrock在x点的梯度grad
% 分别对x1和x2求导
%   此处显示详细说明

grad = [200*(x(2)-x(1).^2)*(-2*x(1)) + 2*(x(1)-1), 200*(x(2)-x(1).^2)];

end

