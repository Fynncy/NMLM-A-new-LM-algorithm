clc;
clear;
close all;

% Sampling
N = 10000; % Sample s i z e
d = 9 ; % Dimension
p = sobolset (d * 2 ) ; % Generating Sobol quasi?random point s e t s
% Rejecting the f i r s t line , taking the l i n e s of 2 to N+1
x0 = p ( 2 :N+1 ,:);
A = x0 ( : , 1 : d ) ; % Matrix A
B = x0 ( : , d+1:d * 2 ) ; % Matrix B
% Matrix A_B^( i )
AB = zeros (N, d , d );
for ii = 1 : d
    AB( : , : , ii ) = A;
    AB( : , ii , ii ) = B( : , ii ) ;
end

YA = myfunc3(A);
YB = myfunc3(B);
YAB = zeros(N,1,d);
for ii = 1:d
    YAB(:,:,ii) = myfunc3(AB(:,:,ii));
end

Var_Y = var(YA) ;
S1 = zeros(1,d);
for ii = 1:d
    sumvar = 0;
    for jj = 1:N
        sumvar = sumvar + YB(jj)*(YAB(jj,1,ii)-YA(jj));
    end
    S1(1,ii) = sumvar/N/Var_Y;
end

ST = zeros(1,d);
for ii = 1:d
    sumvar = 0;
    for jj = 1:N
        sumvar = sumvar + (YA(jj)-YAB(jj,1,ii))^2;
    end
    ST(1,ii) = sumvar/2/N/Var_Y;
end
sum_S1 = sum(S1);
sum_ST = sum(ST);


function y = myfunc(x)

    x(:,1) = x(:,1) * 4 +(-2);
    x(:,2) = x(:,2) * 4 +(-2);
    x(:,3) = x(:,3) * 4 +(-5);
    x(:,4) = x(:,4) * 0.4 +(0.1);
    x(:,5) = x(:,5) * 0.4 +(0.1);
    x(:,6) = x(:,6) * 0.4 +(0.1);
    x(:,7) = x(:,7) * 360 +(-180);
    x(:,8) = x(:,8) * 360 +(-180);
    x(:,9) = x(:,9) * 360 +(-180);
    y = ones(size(x,1),1);
    for jj = 1:size(x,1)
       y(jj,1) = GeneralModel(x(jj,:)); 
    end
end
function [Bs] = GeneralModel( x )
    
    Coil=struct();
    Coil.I=20;       % 线圈电流
    Coil.R=0.4;     % 线圈半径/m
    Coil.f=1000;     % 信号频率
    Coil.Postion=[0 0 0];   %线圈位置
    
    Target=struct();
    Target.Postion=[x(1) x(2) x(3)];  % 目标位置（线圈中心为坐标原点）
    Target.MagPolar=[x(4) x(5) x(6)];  % 目标三轴磁极化率,依次x,y,z
    Target.Theta=x(7)*pi/180;            % 俯仰角        
    Target.Phi=x(8)*pi/180;              % 滚转角 ps:所有角度为0时，代表三轴极化方向与线圈坐标系的x,y,z重合
    Target.Psi=x(9)*pi/180;              % 航向角（垂直放置的物体是没有的）
    
    [Bx, By, Bz] = SecondField(Coil,Target,struct('Postion', Coil.Postion));
%     Bs = sqrt(abs(Bx)^2 + abs(By)^2 + abs(Bz)^2);
    Bs = abs(Bx);
end

function y = myfunc2(x)

    x(:,1) = x(:,1) * 59975000+625000;
    x(:,2) = x(:,2) * 800*10^(-6) +(4*pi*10^(-7));
    x(:,3) = x(:,3) * 0.1 +(0.05);
    x(:,4) = x(:,4) * 0.6 +(0.4);
    x(:,5) = x(:,5) * 9000 +(1000);
%     x(:,2) = x(:,2) * 0.1 +(0.05);
%     x(:,3) = x(:,3) * 0.8 +(0.4);
%     x(:,4) = x(:,4) * 15000 +(500);
    y = ones(size(x,1),1);
    for jj = 1:size(x,1)
       y(jj,1) = PolarizabilityModel2(x(jj,:)); 
    end
end

function beta = PolarizabilityModel( x )

    e = x(1); r = x(2); l = x(3); f = x(4);
    p = l/(2*r); u=4*pi*10^(-7); ur = 1;
    A = 2*pi*r^3*p; %wz19.12.19
    w=2*pi*f;
    %tao=(R^2)*e*u;
    %tao=(R^2)*e*u0/(u);
    tao=(r^2)*e*u/(ur^2); %wz19.12.19
    z1=sqrt(1j*w*tao)-2;
    z2=sqrt(1j*w*tao)+1;
    z0=z1/z2;
    z01=sqrt(1j*w*31*tao)-2;
    z02=sqrt(1j*w*31*tao)+1;
    z00=z01/z02;
    %z(1)=p/2*((1.35-1)+z0);
    %z(2)=p/2*((1.35-1)+z0);
    %z(3)=0.5*((0.3-1)+z00);

    z(1)=0.5*A*((1.05-1)+z0);
    z(2)=0.5*A*((1.05-1)+z0);
    z(3)=p/2*A*((0.9-1)+z00);
    beta = abs(z(2));
end

function beta = PolarizabilityModel2( x )

    e = x(1); u = x(2); r = x(3); l = x(4); f = x(5);
    p = l/(2*r); ur = u/(4*pi*10^(-7)); 
    A = 2*pi*r^3*p; %wz19.12.19
    w=2*pi*f;
    %tao=(R^2)*e*u;
    %tao=(R^2)*e*u0/(u);
    tao=(r^2)*e*u/(ur^2); %wz19.12.19
    z1=sqrt(1j*w*tao)-2;
    z2=sqrt(1j*w*tao)+1;
    z0=z1/z2;
    z01=sqrt(1j*w*31*tao)-2;
    z02=sqrt(1j*w*31*tao)+1;
    z00=z01/z02;
    %z(1)=p/2*((1.35-1)+z0);
    %z(2)=p/2*((1.35-1)+z0);
    %z(3)=0.5*((0.3-1)+z00);

    z(1)=0.5*A*((1.05-1)+z0);
    z(2)=0.5*A*((1.05-1)+z0);
    z(3)=p/2*A*((0.9-1)+z00);
    beta = abs(z(3));
end

function y = myfunc3(x)

    x(:,1) = x(:,1) * 4 +(-2);
    x(:,2) = x(:,2) * 4 +(-2);
    x(:,3) = x(:,3) * 4 +(-5);
    x(:,4) = x(:,4) * 0.02;
    x(:,5) = x(:,5) * 0.02;
    x(:,6) = x(:,6) * 0.02;
    x(:,7) = x(:,7) * 0.02;
    x(:,8) = x(:,8) * 0.02;
    x(:,9) = x(:,9) * 0.02;
    y = ones(size(x,1),1);
    for jj = 1:size(x,1)
       y(jj,1) = LinearModel(x(jj,:)); 
    end
end

 %修改
 
function y = myfunc4(x)

    x(:,1) = x(:,1) * 4 +(-2);
    x(:,2) = x(:,2) * 4 +(-2);
    x(:,3) = x(:,3) * 4 +(-5);
    x(:,4) = x(:,4) * 0.02;
    x(:,5) = x(:,5) * 0.02;
    x(:,6) = x(:,6) * 0.02;
    x(:,7) = x(:,7) * 0.02;
    x(:,8) = x(:,8) * 0.02;
    x(:,9) = x(:,9) * 0.02;
    y = ones(size(x,1),1);
    for jj = 1:size(x,1)
       y(jj,1) = LinearModel(x(jj,:)); 
    end
end

function [Bs] = LinearModel( x )
    
    Coil=struct();
    Coil.I=20;       % 线圈电流
    Coil.R=0.4;     % 线圈半径/m
    Coil.F=1000;     % 信号频率
    Position=[0 0 0];   %线圈位置
    
%     Target=struct();
%     Target.Postion=[x(1) x(2) x(3)];  % 目标位置（线圈中心为坐标原点）
%     Target.MagPolar=[x(4) x(5) x(6)];  % 目标三轴磁极化率,依次x,y,z
%     Target.Theta=x(7)*pi/180;            % 俯仰角        
%     Target.Phi=x(8)*pi/180;              % 滚转角 ps:所有角度为0时，代表三轴极化方向与线圈坐标系的x,y,z重合
%     Target.Psi=x(9)*pi/180;              % 航向角（垂直放置的物体是没有的）
    
    B = LinearizationSecondField(Coil,Position,x);
%     Bs = sqrt(abs(Bx)^2 + abs(By)^2 + abs(Bz)^2);
    Bs = abs(B(3));
end