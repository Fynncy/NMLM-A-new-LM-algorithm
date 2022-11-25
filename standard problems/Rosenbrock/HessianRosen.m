function [hessian] = HessianRosen(x)
%HessianRosen 给定包含两个元素的一位向量x，返回Rosenbrock在x点的海森矩阵hessian
%   此处显示详细说明

hessian = [1200*x(1).^2 - 400*x(2) + 2, -400*x(1); 
            -400*x(1), 200];

end

