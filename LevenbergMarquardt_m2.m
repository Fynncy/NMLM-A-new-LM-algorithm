function [i, xk] = LevenbergMarquardt_m2(f, grad, x0, iterations, tol, m_pk)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

xlog = [];
ylog = [];
xdlog = [];
tol2 = ones(1,length(x0))*tol;
xk = x0;
lamda = 0.01;
updateJ = 1;
muk = 1e-04;
old_gnorm = 0;
iter_num = 0;
filename1 = 'D:\博士项目研究\目标探测\论文\优化\反演\LM算法\matlab_result\final_result\result_grad.csv';
fid1 = fopen(filename1, 'at+');
fprintf(fid1,['%s,'],'mLM2');
for i = 1:iterations
    if updateJ == 1
        xlog = [xlog xk'];
        J = GetResidualGrad_1(m_pk,xk);
        d = f(xk);
        ylog = [ylog d];
        H=J'*J;
        if i==1
            e=d;
        end
    end 
    H_lm=H+(lamda*eye(9,9));
    g = grad(xk);
    old_gnorm = norm(g);
    pk = -H_lm\g';
    xk1 = xk+pk';
    e_lm = f(xk1);
    if e_lm<e
        xk = xk1;
        e = e_lm;
        updateJ=1;
        func_xk =  GetResidualGrad_1(m_pk,xk);
        func_xk1 =  GetResidualGrad_1(m_pk,xk1);
        ro = (func_xk'* func_xk - func_xk1'* func_xk1) / (pk' * (lamda * pk - g));
        rk = norm(ro);
        JF = norm(func_xk)*d;
        if rk > 0.0001
            if JF < (0.25 / muk)  % gnorm→JF→gnorm
                muk = muk * 4;
            elseif JF <= (0.75 / muk)
                muk = muk;
            else
                muk = max(muk * 0.25, 1e-08);
            end
        else
            muk = muk * 2;
        end
        lamda = muk * norm(func_xk)' * (d^2);
    else
        updateJ=0;
        lamda=lamda*10;
    end
    g_xk1 = grad(xk1);
    gnorm = norm(g_xk1);
    xdlog = [xdlog norm(xk' - xlog(:,end))];
%     fprintf("method=%s, iter=%d\n","m-LM2",i);
%         fprintf("method=%s, iter=%d, updateJ=%d, gnorm=%f, lamda=%f\n","m-LM2", i, updateJ, gnorm, lamda);

    fprintf("method=%s, iter=%d, updateJ=%d, gnorm=%f, lamda=%f\n","m-LM2", i, updateJ, gnorm, lamda);
    fprintf(fid1,['%d,'],gnorm);
%     disp(gnorm);
%     disp(old_gnorm);
    if abs(gnorm-old_gnorm) < 1e-05
        iter_num = iter_num + 1;
    else
        iter_num = 0;   
    end
    disp(iter_num);
    if iter_num >= 5
        break;
    end
    

%     
% %     j = grad(x);
% %     H = j' * j;
% %     pk = - H * j';
% %     %alpha = StepLength(f, grad, x, 1.0, pk', c2);
% %     x = x + pk';
%     xdlog = [xdlog norm(xk' - xlog(:,end))];
% %     fprintf("method=%s, iter=%d\n","m-LM2",i);
    
    %a = hessian(x);
%     fprintf("iter=%d, x=%f,%f, f(x)=%f, direction=%f, Hk=%f,%f\n",i, x(1), x(2), ylog(end), grad(x)*pk, a(1,1), a(2,2));
%     if xdlog(end) < tol
%         break;
%     end
    if abs(grad(xk1)) < tol2
%         disp(grad(xk1));
        break;
    end
end
std=fclose('all'); 
% fprintf("iter=%d, x=%f,%f, f(x)=%f\n",i, xk(1), xk(2), ylog(end));
% disp('LM1');
end

