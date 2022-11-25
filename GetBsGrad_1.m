% function [grad] = GetBsGrad(detector, m_pk, xk, epsilon)
% %UNTITLED 此处显示有关此函数的摘要
% %   此处显示详细说明
%     grad = zeros(length(xk),3);
%     ei = zeros(length(xk),1);
%     for k = 1:length(xk)
%         ei(k) = 1.0;
%         d = epsilon * ei;
%         grad(k,:) = (abs(LinearizationSecondField(detector,m_pk,xk+d'))-abs(LinearizationSecondField(detector,m_pk,xk)))' / d(k);
% %         a = LinearizationSecondField(detector,m_pk,xk+d');
% %         b = LinearizationSecondField(detector,m_pk,xk);
% %         grad(k,:) = (a-b)' / d(k);
%         ei(k) = 0.0;
%     end
% end

function [grad] = GetBsGrad_1(m_pk,xk)

%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    %{
    Bs = 100*Gs*M*Hp
    Bs是三个方向的二次场强度，3*1向量；Gs是只和x,y,z三个变量有关的函数，3*3矩阵；M是只和m11,m12.m13,m22,m23,m33六个变量有关的函数，3*3矩阵；Hp只和x,y,z有关，3*1向量；
    Gs,M,Hp在下面都有定义,写成矩阵相乘的形式是下面的样子。
    |Bsx|       |3*(-x)*(-x)*r5-r3,     3*(-x)*(-y)*r5,    3*(-x)*(-z)*r5;  |   |m11, m12, m13;|   |20*0.16*3*x*z*r5;    |                                          
    |Bsy| = 100*|3*(-x)*(-y)*r5,        3*(-y)*(-y)*r5-r3, 3*(-y)*(-z)*r5;  | * |m12, m22, m23;| * |20*0.16*3*y*z*r5;    |                                        
    |Bsz|       |3*(-x)*(-z)*r5,        3*(-y)*(-z)*r5,    3*(-z)*(-z)*r5-r3|   |m13, m23, m33 |   |20*0.16*(3*z*z*r5-r3)|     

    现在要求Bs的三个方向的函数对x,y,z,m11,m12.m13,m22,m23,m33的偏导，
%}
    grad = zeros(length(xk),3);
    %下面一共9个变量x,y,z,m11,m12.m13,m22,m23,m33
    x = xk(1); y = xk(2); z = xk(3);
    m11 = xk(4); m22 = xk(5); m33 = xk(6); m12 = xk(7); m13 = xk(8); m23 = xk(9);
    %M和x,y,z无关
    M = [m11, m12, m13; 
         m12, m22, m23; 
         m13, m23, m33];
     
    x_dt = x - m_pk(1); y_dt = y - m_pk(2); z_dt = z - m_pk(3);
    x_td = m_pk(1) - x; y_td = m_pk(2) - y; z_td = m_pk(3) - z;
    %下面定义Gs，是和x,y,z有关的
    r = sqrt(x_td^2 + y_td^2 + z_td^2);
    r5 = 1 / r^5;%disp('yes');disp(r5);disp('yes');
    r3 = 1 / r^3;
    Gs = [3*x_td*x_td*r5-r3,     3*x_td*y_td*r5,    3*x_td*z_td*r5; 
          3*x_td*y_td*r5,        3*y_td*y_td*r5-r3, 3*y_td*z_td*r5; 
          3*x_td*z_td*r5,        3*y_td*z_td*r5,    3*z_td*z_td*r5-r3]; 

    %下面定义Hp，只和x,y,z有关
    Hp = (1/4)*[20*0.4^2*3*x_dt*z_dt*r5; 
                20*0.4^2*3*y_dt*z_dt*r5; 
                20*0.4^2*(3*z_dt*z_dt*r5-r3)];

    %下面对m11,m22,m33,m12,m23,m33求导，因为只有M中有m11,m22,m33,m12,m23,m33,所以对它们求导的时候Gs和Hp看做常数，后面直接乘就好了
    M1_grad = [1 0 0; 0 0 0; 0 0 0];
    M2_grad = [0 0 0; 0 1 0; 0 0 0];
    M3_grad = [0 0 0; 0 0 0; 0 0 1];
    M4_grad = [0 1 0; 1 0 0; 0 0 0];
    M5_grad = [0 0 1; 0 0 0; 1 0 0];
    M6_grad = [0 0 0; 0 0 1; 0 1 0];
    
    %计算Bs对m11,m22,m33,m12,m23,m33的偏导
    grad(4,:) = 100 * Gs * M1_grad * Hp;
    grad(5,:) = 100 * Gs * M2_grad * Hp;
    grad(6,:) = 100 * Gs * M3_grad * Hp;
    grad(7,:) = 100 * Gs * M4_grad * Hp;
    grad(8,:) = 100 * Gs * M5_grad * Hp;
    grad(9,:) = 100 * Gs * M6_grad * Hp;
    
    %下面是重点，对x,y,z求偏导，先算r对x的偏导，再算r5对x的偏导，再算r3对x的偏导，再算Gs对x的偏导
    x_grad_r = -(x_td) * (x_td^2 + y_td^2 + z_td^2)^(-1/2); %disp(x_grad_r);
    x_grad_r5 = (-5 * r^4 * x_grad_r) / r^10;%disp(x_grad_r5);
    x_grad_r3 = (-3 * r^2 * x_grad_r) / r^6;%disp(x_grad_r3);
    x_grad_Gs = [3*(2*(-x_td)*r5+x_grad_r5*x_td^2)-x_grad_r3, 3*(y_td)*(-r5+x_grad_r5*(x_td)), 3*(z_td)*(-r5+x_grad_r5*(x_td));
                 3*(y_td)*(-r5+x_grad_r5*(x_td)),            3*y_td^2*x_grad_r5-x_grad_r3,    3*(y_td)*(z_td)*x_grad_r5;
                 3*(z_td)*(-r5+x_grad_r5*(x_td)),            3*(y_td)*(z_td)*x_grad_r5,       3*z_td^2*x_grad_r5-x_grad_r3];

    %同样这是对y的
    y_grad_r = -(y_td) * (x_td^2 + y_td^2 + z_td^2)^(-1/2);
    y_grad_r5 = (-5 * r^4 * y_grad_r) / r^10;
    y_grad_r3 = (-3 * r^2 * y_grad_r) / r^6;
    y_grad_Gs = [3*x_td^2*y_grad_r5-y_grad_r3,    3*(x_td)*(-r5+y_grad_r5*(y_td)),            3*(x_td)*(z_td)*y_grad_r5;
                 3*(x_td)*(-r5+y_grad_r5*(y_td)), 3*(2*(-y_td)*r5+y_grad_r5*y_td^2)-y_grad_r3, 3*(z_td)*(-r5+y_grad_r5*(y_td));
                 3*(x_td)*(z_td)*y_grad_r5,       3*(z_td)*(-r5+y_grad_r5*(y_td)),            3*z_td^2*y_grad_r5-y_grad_r3];

    %对z的
    z_grad_r = -(z_td) * (x_td^2 + y_td^2 + z_td^2)^(-1/2);
    z_grad_r5 = (-5 * r^4 * z_grad_r) / r^10;
    z_grad_r3 = (-3 * r^2 * z_grad_r) / r^6;
    z_grad_Gs = [3*x_td^2*z_grad_r5-z_grad_r3,    3*(x_td)*(y_td)*z_grad_r5,       3*(x_td)*(-r5+z_grad_r5*(z_td));
                 3*(x_td)*(y_td)*z_grad_r5,       3*y_td^2*z_grad_r5-z_grad_r3,    3*(y_td)*(-r5+z_grad_r5*(z_td));
                 3*(x_td)*(-r5+z_grad_r5*(z_td)), 3*(y_td)*(-r5+z_grad_r5*(z_td)), 3*(2*(-z_td)*r5+z_grad_r5*z_td^2)-z_grad_r3];

    %这是Hp对x,y,z的偏导
    x_grad_Hp = (1/4)*20*0.4^2* [3*z_dt*(r5+x_grad_r5*x_dt); 3*y_dt*z_dt*x_grad_r5;      3*z_dt^2*x_grad_r5-x_grad_r3];
    y_grad_Hp = (1/4)*20*0.4^2* [3*x_dt*z_dt*y_grad_r5;      3*z_dt*(r5+y_grad_r5*y_dt); 3*z_dt^2*y_grad_r5-y_grad_r3];
    z_grad_Hp = (1/4)*20*0.4^2* [3*x_dt*(r5+z_grad_r5*z_dt); 3*y_dt*(r5+z_grad_r5*z_dt); 3*(2*z_dt*r5+z_grad_r5*z_dt^2)-z_grad_r3];
    
    %因为Gs和Hp都和x,y,z有关，所以要用到乘法求导法则
    grad(1,:) = 100 * (x_grad_Gs*M*Hp + Gs*M*x_grad_Hp);
    grad(2,:) = 100 * (y_grad_Gs*M*Hp + Gs*M*y_grad_Hp);
    grad(3,:) = 100 * (z_grad_Gs*M*Hp + Gs*M*z_grad_Hp);
    
    B_temp = LinearizationSecondField_1(m_pk,xk);
    
    if B_temp(1) < 0
        grad(:,1) = -grad(:,1);
    end
    if B_temp(2) < 0
        grad(:,2) = -grad(:,2);
    end
    if B_temp(3) < 0
        grad(:,3) = -grad(:,3);
    end
end