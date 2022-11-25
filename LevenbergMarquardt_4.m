function [i, xk] = LevenbergMarquardt_4(f, grad, x0, iterations, tol, m_pk)
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
fprintf(fid1,['%s,'],'LM4');
for i = 1:iterations
%     if updateJ == 1
    xlog = [xlog xk'];
    J = GetResidualGrad_1(m_pk,xk);
    d = f(xk);
    ylog = [ylog d];
    H=J'*J;
    if i==1
        e=d;
    end
%     end 
    H_lm=H+(lamda*eye(9,9));
    g = grad(xk);
    pk = -H_lm\g';
    xk1 = xk+pk';
    e_lm = f(xk1);
    if e_lm<e
        e = e_lm;
%         lamda=lamda/10;
%         xk = xk1;%         
%         updateJ=1;
% %     else
% %         updateJ=0;
% %         lamda=lamda*10;
    end
    func_xk =  GetResidualGrad_1(m_pk,xk);
    func_xk1 =  GetResidualGrad_1(m_pk,xk1);
    ro = (func_xk'* func_xk - func_xk1'* func_xk1) / (pk' * (lamda * pk - g));
    rk = norm(ro);
%     disp(rk);
    if rk > 0.0001
        xk = xk1;
    end
    if rk < 0.25
        lamda = lamda * 2;
    elseif rk <= 0.75
        lamda = lamda;
    else
        lamda = lamda * 0.5;
    end
    
    if abs(grad(xk1)) < tol2
%         disp(grad(xk1));
        break;
    end    
    gnorm = norm(g);
    xdlog = [xdlog norm(xk' - xlog(:,end))];
%     fprintf("method=%s, iter=%d\n","LM4",i);
    fprintf("method=%s, iter=%d, updateJ=%d, gnorm=%f, lamda=%f\n","LM4", i, updateJ, gnorm, lamda);
    fprintf(fid1,['%d,'],gnorm);

end
fprintf(fid1,['%s\n,'],'');
std=fclose('all'); 
% fprintf("iter=%d, x=%f,%f, f(x)=%f\n",i, xk(1), xk(2), ylog(end));
% disp('LM4');
end

