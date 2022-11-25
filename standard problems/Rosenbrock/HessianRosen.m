function [hessian] = HessianRosen(x)
%HessianRosen ������������Ԫ�ص�һλ����x������Rosenbrock��x��ĺ�ɭ����hessian
%   �˴���ʾ��ϸ˵��

hessian = [1200*x(1).^2 - 400*x(2) + 2, -400*x(1); 
            -400*x(1), 200];

end

