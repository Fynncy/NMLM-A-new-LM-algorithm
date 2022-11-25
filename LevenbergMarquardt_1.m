function [i, xk] = LevenbergMarquardt_1(f, grad, x0, iterations, tol, m_pk)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

xlog = [];
ylog = [];
xdlog = [];
tol2 = ones(1,length(x0))*tol;
xk = x0;
lamda = 0.01;
alpha = 0;
updateJ = 1;
% wolf求步长
old_fval = f(x0);
old_old_fval = old_fval + norm(grad(x0)) / 2;

filename1 = 'D:\博士项目研究\目标探测\论文\优化\反演\LM算法\matlab_result\final_result\result_lamda.csv';
fid1 = fopen(filename1, 'at+');
fprintf(fid1,['%s,'],'LM1');

function [alpha,xkp1,pkp1,gfkp1,gnorm] = polak_ribiere_powell_step(alpha)
%     tic
    xkp1 = xk + alpha * pk;
    gfkp1 = grad(xkp1);
    yk = gfkp1 - gk;
    beta_k = max(0, (yk*gfkp1')/(gk*gk'));
    pkp1 = -gfkp1 + beta_k * pk;
    gnorm = max(abs(gfkp1));
%     toc
end
function n = descent_condition(alpha, xkp1, fp1, gfkp1)
    [cs1,cs2,cs3,cs4,cs5] = polak_ribiere_powell_step(alpha);
    alphak=cs1;xk1=cs2;pk1=cs3;gk1=cs4;gnorm=cs5;
    if gnorm <= 1e-5
        n = true;
    else
        n = (pk1*gk1') <= -0.01*(gk1*gk1');
    end
end

for i = 1:iterations
    if updateJ == 1
        xlog = [xlog xk'];
%         ylog = [ylog f(xk)];
%         t1 = clock;
        J = GetResidualGrad_1(m_pk,xk);
%         t2 = clock;
%         fprintf("j=%f\n",etime(t2,t1));
        d = f(xk); % d → yk
        ylog = [ylog d];
        H=J'*J;
        if i==1
            e=d;
        end
    end 
    
    old_fval_back = old_fval; %wz-19.12.16
    old_old_fval_back = old_old_fval;
    
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
%         lamda=lamda/10;
        xk = xk1;
        e = e_lm;
%         disp(e);
        updateJ=1;
        lamda=(norm(d))^2;
    else
%         % 需要修改的部分
% %         [alpha, old_fval, old_old_fval] = StepLength3(f, grad, xk, pk', grad(xk), old_fval, old_old_fval, 1e-4, 0.4, 1e100, 1e-100, 1e-14);
%         
%         tic
%         [alpha, old_fval, old_old_fval] = StepLength3(f, grad, xk, pk', g, old_fval, old_old_fval, 1e-4, 0.9, 1e100, 1e-100, 1e-14);
%         if(isnan(alpha))
%             [alpha, old_fval, old_old_fval] = StepLength2(f, grad, xk, pk', grad(xk), old_fval_back, old_old_fval_back, 1e-4, 0.9, 1e100, @descent_condition, 1, 20);
%         end
%         toc        
% %         xk = xk + alpha * pk;
        updateJ=0;
        lamda=lamda*10;
%         lamda=(norm(d))^2;
%         fprintf("iter=%d, updateJ=%d, alpha=%f\n",i, updateJ, alpha);
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
    %a = hessian(x);
%     fprintf("iter=%d, x=%f,%f, f(x)=%f, direction=%f, Hk=%f,%f\n",i, x(1), x(2), ylog(end), grad(x)*pk, a(1,1), a(2,2));
%     if xdlog(end) < tol
%         break;
%     end
%     fprintf("iter=%d, lamda=%f,updateJ=%d, alpha=%f\n",i, lamda, updateJ, alpha);
%     fprintf("method=%s, iter=%d\n","LM1",i);
    fprintf("method=%s, iter=%d, updateJ=%d, gnorm=%f, lamda=%f\n","LM1", i, updateJ, gnorm, lamda);
    fprintf(fid1,['%d,'],gnorm);  
end
fprintf(fid1,['%s\n,'],'');
std=fclose('all'); 
% fprintf("iter=%d, x=%f,%f, f(x)=%f\n",i, xk(1), xk(2), ylog(end));
% disp('LM1');
end

