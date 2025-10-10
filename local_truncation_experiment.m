% WORK IN PROGRESS. I think my function calls are wrong and 
% I need to incorporate loglog_fit but not sure how
% MAYBE [p,k] = loglog_fit(h_step, errs, 1)

% Local: the error for a single time step
% log-log fit to estimate order p

%[p,k]: the regressed values for relationship y = k*x^p
%function [p,k] = loglog_fit(x_regression,y_regression,varargin)

function local_truncation_experiment()
    tref = 0.492;
    sol = @solution01;
    rate = @rate_func01;

    % Single step function methods
    forEu = @forward_euler_step
    expMid = @explicit_midpoint_step

    h_step = logspace(-5, 1, 100) % 1e-5 to 1e1 100 points
    figure; hold on;


        errs = zeros(size(h_step));
        for i = 1:length(h_step)
            h = h_step(i);
            XA = sol(tref); % finding exact X at the tref X spot
            [XB, ~] = forEu(rate, tref, XA, h);  % one step from exact state
            X_true = sol(tref + h);
            errs(i) = norm(XB - X_true);        % norm handles vector/scalar

        end
        [p,k] = loglog_fit(h_step(:), errs(:)); % converts to column vectors
        loglog(h_step, errs, 'o-'); hold on;
        disp(p);
        disp(k);
        xlabel('h (step size)');
        ylabel('Error');
        legend('Forward Euler', 'Explicit Midpoint');
% 
        for i = 1:length(h_step)
            h = h_step(i);
            XA = sol(tref); % finding exact X at the tref X spot
            [XB, ~] = expMid(rate, tref, XA, h);  % one step from exact state
            X_true = sol(tref + h);
            errs(i) = norm(XB - X_true);        % norm handles vector/scalar
           
        end
        [p,k] = loglog_fit(h_step(:), errs(:));
        loglog(h_step, errs, 'o-');
        disp(p);
        disp(k);
        xlabel('h (step size)');
        ylabel('Error');
        legend('Forward Euler', 'Explicit Midpoint');

        % is y = k*x^p the same as e = O*h^p? 
        % [p,k] = loglog_fit(h_step, errs, 1)

end