% function [fval] = ObjectFun3(Alpha,res)
% %UNTITLED4 �˴���ʾ�йش˺�����ժҪ
% %   �˴���ʾ��ϸ˵��
%     if(size(Alpha,1)>1)
%         Alpha = Alpha';
%     end
%     v_Residual = res(Alpha);
%     fval = sum(v_Residual.^2)/2;
% end

function [fval] = ObjectFun3_1(Alpha,m_pk,m_data)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    if(size(Alpha,1)>1)
        Alpha = Alpha';
    end
    v_Residual = GetResidualVector3_1(Alpha,m_pk,m_data);
    fval = sum(v_Residual.^2)/2;
end
