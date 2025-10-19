% CURRENT

% Local: the error for a single time step
% log-log fit to estimate order p

%[p,k]: the regressed values for relationship y = k*x^p
%function [p,k] = loglog_fit(x_regression,y_regression,varargin)

function local_truncation_experiment1()
    tref = 0.492;
    sol = @solution01;
    rate = @rate_func01;

    % Single step function methods
    forEu = @forward_euler_step
    expMid = @explicit_midpoint_step

    backEu = @backward_euler_step
    impMid = @implicit_midpoint_step


    h_step = logspace(-5, 1, 100); % 1e-5 to 1e1 100 points
    
    analytical_diff  = zeros(size(h_step));
    errs_FE = zeros(size(h_step));
    errs_Mid = zeros(size(h_step));

    XA = sol(tref); % finding exact X at the tref X spot

    for r = 1:length(h_step)
        h = h_step(r);

        
        [XB_FE, ~] = forEu(rate, tref, XA, h);  % one step from exact state
        [XB_Mid, ~] = expMid(rate, tref, XA, h);  % one step from exact state
        [XB_backEu, ~] = backEu(rate, tref, XA, h);  % one step from exact state
        [XB_impMid, ~] = impMid(rate, tref, XA, h);  % one step from exact state
        X_true = sol(tref + h);

        analytical_diff(r) = norm(X_true - XA);
        errs_FE(r) = norm(XB_FE - X_true);
        errs_Mid(r) = norm(XB_Mid - X_true);
        errs_backEu(r) = norm(XB_backEu - X_true);
        errs_impMid(r) = norm(XB_impMid - X_true);

        % errs(i) = norm(XB - X_true);        % norm handles vector/scalar

    end
    % [p,k] = loglog_fit(h, errs(i), 1);
    % plot(errs(i), h);
    
    figure(1);
%     loglog(h_step,analytical_diff,'ko','MarkerFaceColor','k','MarkerSize',1);
%     hold on
    loglog(h_step,errs_FE,'ro','MarkerFaceColor','r','MarkerSize',1); hold on;
    loglog(h_step,errs_Mid,'bo','MarkerFaceColor','b','MarkerSize',1);

    loglog(h_step,errs_backEu,'go','MarkerFaceColor','g','MarkerSize',1);
    loglog(h_step,errs_impMid,'mo','MarkerFaceColor','m','MarkerSize',1);


    filter_params = struct();

    filter_params.max_xval = 1;

%     [p1,k1] = loglog_fit(h_step, analytical_diff, filter_params);
    [p2,k2] = loglog_fit(h_step, errs_FE, filter_params);
    [p3,k3] = loglog_fit(h_step, errs_Mid, filter_params);
    
    [p4,k4] = loglog_fit(h_step, errs_backEu, filter_params);
    [p5,k5] = loglog_fit(h_step, errs_impMid, filter_params);
% 
%     loglog(h_step,k1*h_step.^p1,'g','LineWidth',1)
%     loglog(h_step,k2*h_step.^p2,'r','LineWidth',1)
%     loglog(h_step,k3*h_step.^p3,'b','LineWidth',1)
% 
%     loglog(h_step,k4*h_step.^p4,'y','LineWidth',1)
%     loglog(h_step,k5*h_step.^p5,'m','LineWidth',1)


    legend(sprintf('Forward Euler (p = %.4f)', p2), sprintf('Explict Midpoint (p = %.4f)', p3),  sprintf('Backward Euler (p = %.4f)', p4),  sprintf('Implicit Midpoint (p = %.4f)', p5));
    xlabel('h step')
    ylabel('Error')

    title('Local Truncation Error for All Methods to Solution 1')

    disp(p1)
    disp(p2)
    disp(p3)
    disp(p4)
    disp(p5)

    xlabel('')
    ylabel('Error')

    % is y = k*x^p the same as e = O*h^p? 
    % [p,k] = loglog_fit(h_step, errs, 1)
 
end