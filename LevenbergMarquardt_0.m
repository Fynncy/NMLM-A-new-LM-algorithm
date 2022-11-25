function [i, xk] = LevenbergMarquardt_0(f, grad, x0, iterations, tol, m_pk)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

xlog = [];
ylog = [];
xdlog = [];
tol2 = ones(1,length(x0))*tol;
xk = x0;
lamda = 0.01;
updateJ = 1;
filename1 = 'D:\博士项目研究\目标探测\论文\优化\反演\LM算法\matlab_result\final_result\result_lamda.csv';
fid1 = fopen(filename1, 'at+');
fprintf(fid1,['%s,'],'LM0');
for i = 1:iterations
    if updateJ == 1
        xlog = [xlog xk'];
%         ylog = [ylog f(xk)];
%         t1 = clock;
        J = GetResidualGrad_1(m_pk,xk);
%         t2 = clock;
%         fprintf("j=%f\n",etime(t2,t1));
        d = f(xk);
        ylog = [ylog d];
        H=J'*J;
        if i==1
            e=d;
        end
    end 
    H_lm=H+(lamda*eye(9,9));
%     t3 = clock;
    g = grad(xk);
%     t4 = clock;
%     fprintf("g=%f\n",etime(t4,t3));
%     t5 = clock;
    pk = -H_lm\g';
%     t6 = clock;
%     fprintf("pk=%f\n",etime(t6,t5));
    xk1 = xk+pk';
    e_lm = f(xk1);
    if e_lm<e
        lamda=lamda/10;
        xk = xk1;
        e = e_lm;
%         disp(e);
        updateJ=1;
    else
        updateJ=0;
        lamda=lamda*10;
    end
    if abs(grad(xk1)) < tol2
%         disp(grad(xk1));
        break;
    end
    gnorm = norm(g);
    
%     j = grad(x);
%     H = j' * j;
%     pk = - H * j';
%     %alpha = StepLength(f, grad, x, 1.0, pk', c2);
%     x = x + pk';
    xdlog = [xdlog norm(xk' - xlog(:,end))];
%     fprintf("method=%s, iter=%d\n","LM0",i);
%     fprintf("iter=%d, updateJ=%d, gnorm=%f, lamda=%f\n",i, updateJ, gnorm, lamda);
    fprintf("method=%s, iter=%d, updateJ=%d, gnorm=%f, lamda=%f\n","LM0", i, updateJ, gnorm, lamda); 
%     fprintf(fid1,['%d,'],lamda); 
    fprintf(fid1,['%d,'],gnorm); 
    %a = hessian(x);
%     fprintf("iter=%d, x=%f,%f, f(x)=%f, direction=%f, Hk=%f,%f\n",i, x(1), x(2), ylog(end), grad(x)*pk, a(1,1), a(2,2));
%     if xdlog(end) < tol
%         break;
%     end
end
fprintf(fid1,['%s\n,'],'');
std=fclose('all'); 
% fprintf("iter=%d, x=%f,%f, f(x)=%f\n",i, xk(1), xk(2), ylog(end));
% disp('LM1');
end

