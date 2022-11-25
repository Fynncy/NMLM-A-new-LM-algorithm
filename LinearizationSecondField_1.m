function [B] = LinearizationSecondField_1(Pk,Alpha)
%������Ȧ����������ģ�Ͳ���������γ� ���������۲ο���֪���������-2.3�ڣ�
%   ���: 
%       B��3*1���������γ�
%   ����: 
%       detector��������������һ���ṹ�壬������Ȧ�������뾶���Լ�����Ƶ��
%       pk,������λ��
%       Alpha��1*9����������Ϊ��λ�ò���x,y,z��6����Ԫ��M11,M22,M33,M12,M13,M23 
  v_r=Alpha(1:3);
  v_M=Alpha(4:9);
  v_r=v_r(:);
  v_M=v_M(:);
  Pk=Pk(:);
  [~,~,~,Gk]=HFieldModel([v_r v_r],Pk(1),Pk(2),Pk(3));      % �����ʽ����Gk 
%   Bp=FirstField(detector.I,detector.R,detector.F,v_r-Pk);   % ����̽��Ŀ�괦������Bp
  Bp=FirstField_wz_matrix(20,0.4,(v_r-Pk)');   % ����̽��Ŀ�괦������Bp %wz19.12.19
  Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...                           % ����Wk
      0 Bp(2) 0 Bp(1) 0 Bp(3);...
      0 0 Bp(3) 0 Bp(1) Bp(2)];
   B=Gk*Wk*v_M;
end

