function [grad] = JacobianRosen(x)
%JacobianRosen ������������Ԫ�ص�һλ����x������Rosenbrock��x����ݶ�grad
% �ֱ��x1��x2��
%   �˴���ʾ��ϸ˵��

grad = [200*(x(2)-x(1).^2)*(-2*x(1)) + 2*(x(1)-1), 200*(x(2)-x(1).^2)];

end

