function [y] = Rosenbrock(x)
%Rosenbrock ������������Ԫ�ص�һλ����x������Rosenbrock��x��ĺ���ֵy
%   �˴���ʾ��ϸ˵��

y = 100 .* (x(2) - x(1).^2).^2 + (1 - x(1)).^2;

end

