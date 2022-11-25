function [B] = LinearizationSecondField_1(Pk,Alpha)
%给定线圈参数，根据模型参数计算二次场 （具体理论参考刘知洋毕设论文-2.3节）
%   输出: 
%       B，3*1向量，二次场
%   输入: 
%       detector，发射器参数，一个结构体，包括线圈电流、半径和以及发射频率
%       pk,发射器位置
%       Alpha，1*9向量，依次为三位置参数x,y,z和6张量元素M11,M22,M33,M12,M13,M23 
  v_r=Alpha(1:3);
  v_M=Alpha(4:9);
  v_r=v_r(:);
  v_M=v_M(:);
  Pk=Pk(:);
  [~,~,~,Gk]=HFieldModel([v_r v_r],Pk(1),Pk(2),Pk(3));      % 计算格式函数Gk 
%   Bp=FirstField(detector.I,detector.R,detector.F,v_r-Pk);   % 计算探测目标处的主场Bp
  Bp=FirstField_wz_matrix(20,0.4,(v_r-Pk)');   % 计算探测目标处的主场Bp %wz19.12.19
  Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...                           % 计算Wk
      0 Bp(2) 0 Bp(1) 0 Bp(3);...
      0 0 Bp(3) 0 Bp(1) Bp(2)];
   B=Gk*Wk*v_M;
end

