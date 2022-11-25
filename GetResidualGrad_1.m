% function [grad] = GetResidualGrad(xk,f,epsilon)
% %UNTITLED 此处显示有关此函数的摘要
% %   此处显示详细说明
%     grad = [];
%     ei = zeros(length(xk),1);
%     fxk = f(xk);
%     for k = 1:length(xk)
%         ei(k) = 1.0;
%         d = epsilon * ei;
%         grad(:,k) = (f(xk+d')-fxk) / d(k);
%         ei(k) = 0.0;
%     end
% end

function [grad] = GetResidualGrad_1(m_pk,xk)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    grad = zeros(3*size(m_pk,1),9);
    for k = 1:size(m_pk,1)
%         BsGrad = GetBsGrad(detector,m_pk(k,:),xk,1.5e-8);
        BsGrad = GetBsGrad_1(m_pk(k,:),xk);
        grad(3*k-2,:) = BsGrad(:,1)';
        grad(3*k-1,:) = BsGrad(:,2)';
        grad(3*k,:) = BsGrad(:,3)';
    end
end


