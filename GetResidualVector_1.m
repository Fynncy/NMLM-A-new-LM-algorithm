function [v_Residual] = GetResidualVector_1(Alpha,m_pk,m_data)
%   残差函数
%   输出: 
%       v_Residual，残差向量,mx1向量
%   输入: 
%       detector，发射器参数
%       Alpha，1*9向量，依次为x,y,z,M11,M22,M33,M12,M13,M23
%       m_pk,m个发射器位置
%       m_data,m组观测数据
v_Residual=zeros(3*size(m_pk,1),1);      % 残差向量，长度需同所用数据组数一致。
for i=1:size(m_pk,1)
    B=LinearizationSecondField_1(m_pk(i,:),Alpha);
%     v_Residual(i)=norm(m_data(i,1:3),2)-norm(B(1:3)',2);
%% 仅用z轴观测值来构建残差函数
%       v_Residual(i)=m_data(i,3)-B(3);
      
%% 同时用x,y,z轴观测值来构建残差函数
    v_Residual(3*i-2)=m_data(i,1)-B(1);
    v_Residual(3*i-1)=m_data(i,2)-B(2);
    v_Residual(3*i)=m_data(i,3)-B(3);
end
%   fprintf('alpha = %d %d %d %d %d %d %d %d %d\n',Alpha(1),Alpha(2),Alpha(3),Alpha(4),Alpha(5),Alpha(6),Alpha(7),Alpha(8),Alpha(9));
%   disp(v_Residual(1:9));
end

