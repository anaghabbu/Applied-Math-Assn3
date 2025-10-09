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

    for i = 1:length(forEu)
        errs = zeros(size(h_step));
        for r = 1:length(h_step)
            h = h_step(i)
            XA = sol(tref) % finding exact X at the tref X spot
            [XB, ~] = forEu(rate, tref, XA, h);  % one step from exact state
            [XB, ~] = expMid(rate, tref, XA, h);  % one step from exact state
            X_true = sol(tref + h);
            errs(i) = norm(XB - X_true);        % norm handles vector/scalar
            [p,k] = loglog_fit(h, errs(i), 1);
            plot(errs(i), h);
            
        end
        % is y = k*x^p the same as e = O*h^p? 
        % [p,k] = loglog_fit(h_step, errs, 1)
    end   
end