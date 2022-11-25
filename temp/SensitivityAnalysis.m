clc;
clear;
close all;

% Sampling
N = 10000; % Sample s i z e
d = 4 ; % Dimension
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

YA = myfunc(A);
YB = myfunc(B);
YAB = zeros(N,1,d);
for ii = 1:d
    YAB(:,:,ii) = myfunc(AB(:,:,ii));
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


disp(S1);
disp(sum_S1);
disp(ST);
disp(sum_ST);
% fprintf(sum_S1);
% fprintf(sum_ST);
disp('finsh')



function y = myfunc(x)
    x(:,1) = x(:,1) * 9 +(-5);
    x(:,2) = x(:,2) * 5 +(-2);
    x(:,3) = x(:,3) * 4 +(-2);
    x(:,4) = x(:,4)* 4 +(-2);
%     x(:,5) = x(:,5) * 0.02;
%     x(:,6) = x(:,6) * 0.02;
%     x(:,7) = x(:,7) * 0.02;
%     x(:,8) = x(:,8) * 0.02;
%     x(:,9) = x(:,9) * 0.02;
    y = ones(size(x,1),1);
    for jj = 1:size(x,1)
       y(jj,1) = Powellf(x(jj,:)); 
    end
end

function y = Powellf(x)
%Rosenbrock 给定包含两个元素的一位向量x，返回Rosenbrock在x点的函数值y
%   此处显示详细说明

% y = 100 .* (x(2) - x(1).^2).^2 + (1 - x(1)).^2;
y = (x(1) + 10.* x(2)).^2 + 5 .* (x(3) - x(4)).^2 + ...
    (x(2) - 2 .* x(3)).^4 + 10 .* (x(1) - x(4)).^4;
% 最小值(0,0,0,0)
end
