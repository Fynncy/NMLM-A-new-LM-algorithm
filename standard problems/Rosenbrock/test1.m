clear;
clc;
close all;
format long

global coordX coordY iter
iter = 0;

% graph initialization
figure(1)
figure(2)

% ==================== Select a starting point ====================
startingpoint = 4; % Choose a starting point between 1-4 by changing the starting point variable
X0 = [2 0.5; 1 -1.5; -1 -1.5; -1 0.5]; % starting points
x0 = [X0(startingpoint,1) X0(startingpoint,2)]; % starting point vector value

figure(1) %  2D contour plot of a function with marked trajectories
hold on
axis tight
[X, Y] = meshgrid(-4:0.001:4,-4:0.001:4);
Z = plotRosenbrock(X, Y); % compute the contour value of the Rosenbrock function
[M, c] = contourf(X,Y,log(Z),'ShowText','on');
c.Fill = 0; % exclusion of colors
c.LineWidth = 0.33; % reduction of line size for better visibility

% If there is no graphs folder, create one.
if ~exist("./graphs", 'dir')
       mkdir("./graphs")
end

% ==================== Calling individual methods ====================
% select a method by changing the value of the choice variable, where:
% 1 - means the Quasi-Newton method
% 2 - means the Region of Trust method without the given Hessian
% 3 - means the Region of Trust method with the given Hessian
% 4 - means the Nelder-Mead method

choice = 4;

switch choice
    case 1
        [x, value] = optimQuasiNewton(x0, startingpoint);
    case 2
        [x, value] = optimTrustRegion(x0, startingpoint);
    case 3
        [x, value] = optimTrustRegionHessian(x0, startingpoint);
    case 4
        [x, value] = optimSimplex(x0, startingpoint);
end


% !!! before running the script, you must install: Optimization Toolbox
% source: https://www.google.com/search?q=Optimization+Toolbox&sourceid=chrome&ie=UTF-8
function [x, value] = optimQuasiNewton(x0, startingpoint) % the function responsible for the optimization of the unconsidered function by the Quasi-Newton method
    
    global coordX coordY iter
    % setting the appropriate settings for the solver
    options = optimoptions(@fminunc,'OutputFcn',@outputQuasiNewton,'Display', ...
            'iter-detailed','Algorithm','quasi-newton','MaxIterations', 1000, 'StepTolerance', ...
            1e-12, 'FunctionTolerance', 1e-12);
    [x, value] = fminunc(@rosenbrock, x0, options); % solver invocation

    function stop = outputQuasiNewton(x, optimValues, state)
        stop = false; % setting the stop variable responsible for the operation of the function
        if isequal(state,'init') % preparation of charts for marking subsequent iterations
            figure(1)
            hold on
            title('Quasi-Newton method')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("x value")
            ylabel("y value")

            figure(2)
            hold on
            title('Quasi-Newton method (objective function values)');
            xlabel("number of iterations [n]")
            ylabel("function value [log(y)]")
            set(gca, 'YScale', 'log') % displaying the Y axis on a logarithmic scale

        elseif isequal(state,'iter') % updating the graph in subsequent iterations
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % support for the final iteration of a function
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = char("./graphs/QuasiNewton1_" + startingpoint + ".png");
            saveas(gcf,path);
            fprintf('Grapgh no. 1 has been saved\n')
            hold off
            figure(2)
            path = char("./graphs/QuasiNewton2_" + startingpoint + ".png");
            saveas(gcf,path);
            fprintf('Grapgh no. 2 has been saved\n')
            hold off 
    
        end
    end
end

% !!! before running the script, you must install: Optimization Toolbox
% source: https://www.google.com/search?q=Optimization+Toolbox&sourceid=chrome&ie=UTF-8
function [x, value] = optimTrustRegion(x0, startingpoint) % the function responsible for the optimization of the unconsidered function by the Trust-Region method without the given Hessian
    
    global coordX coordY iter
    % setting the appropriate settings for the solver
    options = optimoptions(@fminunc,'OutputFcn',@outputTrustRegion,'Display', ...
            'iter-detailed','Algorithm','trust-region', 'SpecifyObjectiveGradient',true, ...
            'MaxIterations', 1000, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12);
    [x, value] = fminunc(@rosenbrockwithgrad, x0, options); % solver invocation
    
    function stop = outputTrustRegion(x, optimValues, state)
        stop = false; % setting the stop variable responsible for the operation of the function
        if isequal(state,'init') % preparation of charts for marking subsequent iterations
            figure(1)
            hold on
            title('Trust-Region method without the given Hessian')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("x value")
            ylabel("y value")

            figure(2)
            hold on
            title('Trust-Region method without the given Hessian (objective function values)');
            xlabel("number of iterations [n]")
            ylabel("function value [log(y)]")
            set(gca, 'YScale', 'log') % displaying the Y axis on a logarithmic scale

        elseif isequal(state,'iter') % updating the graph in subsequent iterations
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % support for the final iteration of a function
            fprintf('done \n')
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = "./graphs/TrustRegion1_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Grapgh no. 1 has been saved\n')
            hold off
            figure(2)
            path = "./graphs/TrustRegion2_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Grapgh no. 2 has been saved\n')
            hold off 
            
        end
    end
end

% !!! before running the script, you must install: Optimization Toolbox
% source: https://www.google.com/search?q=Optimization+Toolbox&sourceid=chrome&ie=UTF-8
function [x, value] = optimTrustRegionHessian(x0, startingpoint) % function responsible for optimizing the considered function using the Trust-Region method with the given Hessiann
    
    global coordX coordY iter
    % setting the appropriate settings for the solver
    options = optimoptions(@fminunc,'OutputFcn',@outputTrustRegionHessian,'Display', ...
            'iter-detailed','Algorithm','trust-region', 'SpecifyObjectiveGradient',true, ...
            'MaxIterations', 1000, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12, 'HessianFcn', 'objective');
    [x, value] = fminunc(@rosenbrockwithhes, x0, options); % solver invocation
    
    function stop = outputTrustRegionHessian(x, optimValues, state)
        stop = false; % setting the stop variable responsible for the operation of the function
        if isequal(state,'init') % preparation of charts for marking subsequent iterations
            figure(1)
            hold on
            title('Trust-Region method with the given Hessiann')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("x value")
            ylabel("y value")

            figure(2)
            hold on
            title('Trust-Region method with the given Hessiann (objective function values)');
            xlabel("number of iterations [n]")
            ylabel("function value [log(y)]")
            set(gca, 'YScale', 'log') % displaying the Y axis on a logarithmic scale

        elseif isequal(state,'iter') % updating the graph in subsequent iterations
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % support for the final iteration of a function
            fprintf('done \n')
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = "./graphs/TrustRegionHessian1_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Grapgh no. 1 has been saved\n')
            hold off
            figure(2)
            path = "./graphs/TrustRegionHessian2_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Grapgh no. 2 has been saved\n')
            hold off 
        end
    end
end

function [x, value] = optimSimplex(x0, startingpoint) % the function responsible for the optimization of the considered function by the Nelder-Mead method
    
    global coordX coordY iter
    
    options = optimset('OutputFcn', @outputSimplex,...
            'MaxFunEvals', 1000, 'TolX', 1e-12, 'TolFun', 1e-12);
    [x, value, exitflag, output] = fminsearch(@rosenbrock, x0, options); % solver invocation
    
    function stop = outputSimplex(x, optimValues, state)
        stop = false; % setting the stop variable responsible for the operation of the function
        if isequal(state,'init') % preparation of charts for marking subsequent iterations
            figure(1)
            hold on
            title('Nelder-Mead method')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("x value")
            ylabel("y value")

            figure(2)
            hold on
            title('Nelder-Mead method (objective function values)');
            xlabel("number of iterations [n]")
            ylabel("function value [log(y)]")
            set(gca, 'YScale', 'log') % displaying the Y axis on a logarithmic scale

        elseif isequal(state,'iter') % updating the graph in subsequent iterations
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % support for the final iteration of a function
            fprintf('done \n')
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = "./graphs/Nelder-Mead1_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Grapgh no. 1 has been saved\n')
            hold off
            figure(2)
            path = "./graphs/Nelder-Mead2_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Grapgh no. 2 has been saved\n')
            hold off 
            
        end
    end
    output
end


% The values of the constants a = 0 and b = -0.5 for the Rosenbrock ("banana") function:
% f(x) = (1-x+a)^2 + 100[y-b-(x-a)^2]^2

function f = plotRosenbrock(x, y) % function required to plot the Rosenbrock function
    f = (1-x).^2 + 100 * (y + 0.5 - x.^2).^2;
end

function f = rosenbrock(x) % proper function used to optimize Rosenbrock function
    f = (1-x(1)).^2 + 100 * (x(2) + 0.5 - x(1).^2).^2;
end

function [f,g] = rosenbrockwithgrad(x) % correct function used to optimize Rosenbrock functions with defined gradient
    f = (1-x(1)).^2 + 100 * (x(2) + 0.5 - x(1).^2).^2;

    if nargout > 1 % if there is more than one argument that the function returns
        g = [400*((-x(2)-0.495)*x(1)+x(1).^3-0.005); % gradient calculated, according to the formula: https://wikimedia.org/api/rest_v1/media/math/render/svg/d632a346cd0677aef80d9fa32f476a5b5bf4dc58
        200*(x(2)-x(1).^2+0.5)];
    end
end

function [f, g, H] = rosenbrockwithhes(x) % function used to optimize a Rosenbrock function with a defined gradient and given hesian
    f = (1 - x(1)).^2 + 100 * (x(2) + 0.5 - x(1).^2).^2;

    if nargout > 1

        g = [400 * ((-x(2) - 0.495) * x(1) + x(1).^3 - 0.005);
        200 * (x(2) - x(1).^2 + 0.5)];

        if nargout > 2

            H = [-400*(-1 * x(1).^2 + x(2) + 0.5) + 800 * x(1).^2 + 2, -400 * x(1);
            -400 * x(1), 200];
        end 
    end
end