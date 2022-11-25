% function [grad] = GetObjGrad3(Alpha,res)
% %UNTITLED3 此处显示有关此函数的摘要
% %   此处显示详细说明
%     if(size(Alpha,1)>1)
%         Alpha = Alpha';
%     end
%     rx = res(Alpha);
%     jx = GetResidualGrad(Alpha,res,1.5e-8);
%     grad = (jx'*rx)';
% end

function [grad] = GetObjGrad3_1(Alpha,m_pk,m_data)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    if(size(Alpha,1)>1)
        Alpha = Alpha';
    end
    rx = GetResidualVector3_1(Alpha,m_pk,m_data);
    jx = GetResidualGrad_1(m_pk,Alpha);
    grad = (jx'*rx)';
end


