function global_truncation_experiment()

    t0 = 0;
    tf = 1;
    sol = @solution01;    
    rate = @rate_func01;

    X0 = sol(t0);  
    X_true = sol(tf);   

    forEu = @forward_euler_fixed_step_integration;
    expMid = @explicit_midpoint_fixed_step_integration;

    h_step = logspace(-5, -0.1, 50); 

    errs_FE  = zeros(size(h_step));
    errs_Mid = zeros(size(h_step));
    nfe_FE   = zeros(size(h_step));
    nfe_Mid  = zeros(size(h_step));

    % Step Size Loop
    for r = 1:length(h_step)
        h = h_step(r);

        % Forward Euler
        [t_FE, X_FE, h_avg_FE, num_evals_FE] = forEu(rate, [t0 tf], X0, h);
        X_FE_final = X_FE(end,:)';
        errs_FE(r) = norm(X_FE_final - X_true);
        nfe_FE(r) = num_evals_FE;

        % Midpoint
        [t_Mid, X_Mid, h_avg_Mid, num_evals_Mid] = expMid(rate, [t0 tf], X0, h);
        X_Mid_final = X_Mid(end,:)';
        errs_Mid(r) = norm(X_Mid_final - X_true);
        nfe_Mid(r)  = num_evals_Mid;
    end

    % Error vs Step Size
    figure(1); clf
    loglog(h_step, errs_FE,  'ro', 'MarkerFaceColor','r', 'MarkerSize',2); hold on
    loglog(h_step, errs_Mid, 'bo', 'MarkerFaceColor','b', 'MarkerSize',2);

    filter_params = struct();
    filter_params.max_xval = 1;

    [p_FE, k_FE]   = loglog_fit(h_step, errs_FE, filter_params);
    [p_Mid, k_Mid] = loglog_fit(h_step, errs_Mid, filter_params);

    loglog(h_step, k_FE*h_step.^p_FE, 'r-', 'LineWidth', 1);
    loglog(h_step, k_Mid*h_step.^p_Mid, 'b-', 'LineWidth', 1);

    legend(sprintf('Forward Euler (p = %.2f)', p_FE), sprintf('Midpoint (p = %.2f)', p_Mid),'Location','northwest');

    xlabel('step size h');
    ylabel('global truncation error');
    title('Explicit Global Truncation Error vs Step Size');
end
