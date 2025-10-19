function explicit_vs_implicit_global()

    t0 = 0;
    tf = 1;
    sol = @solution01;    
    rate = @rate_func01;
    tspan = [0, 10]; % tspan not working in this situation
    h_ref = 0.01;

    X0 = sol(t0);  
    X_true = sol(tf);   


    h_step = logspace(-4, -0.1, 50); 

    % errors

    errs_F  = zeros(size(h_step));
    nfe_F   = zeros(size(h_step));

    errs_B  = zeros(size(h_step));
    nfe_B   = zeros(size(h_step));

    errs_exp_mid = zeros(size(h_step));
    nfe_exp_mid  = zeros(size(h_step));

    errs_imp_mid = zeros(size(h_step));
    nfe_imp_mid  = zeros(size(h_step));

    % Step Size Loop
    for r = 1:length(h_step)
        h = h_step(r);
      
        % Forward Euler
   
        [t_F, X_F, h_avg_F, num_evals_F] = forward_euler_fixed_step_integration(rate, [t0 tf], X0, h);
        
        X_F_final = X_F(end,:)';
        errs_F(r) = norm(X_F_final - X_true);
        nfe_F(r) = num_evals_F;

        % Backward Euler
        
        [t_B, X_B, h_avg_B, num_evals_B] = backward_euler_fixed_step_integration(rate, [t0 tf], X0, h);
 
        X_B_final = X_B(end,:)';
        errs_B(r) = norm(X_B_final - X_true);
        nfe_B(r) = num_evals_B;

        % Explicit Midpoint
    
        [t_exp_mid, X_exp_mid, h_avg_exp_mid, num_evals_exp_mid] = explicit_midpoint_fixed_step_integration(rate, [t0 tf], X0, h);
        X_exp_mid_final = X_exp_mid(end,:)';
        errs_exp_mid(r) = norm(X_exp_mid_final - X_true);
        nfe_exp_mid(r)  = num_evals_exp_mid;

    
        [t_imp_mid , X_imp_mid , h_avg_imp_mid , num_evals_imp_mid ] = implicit_midpoint_fixed_step_integration(rate, [t0 tf], X0, h);
        X_imp_mid_final = X_imp_mid(end,:)';
        errs_imp_mid (r) = norm(X_imp_mid_final - X_true);
        nfe_imp_mid(r)  = num_evals_imp_mid ;
    end

    % Error vs Step Size
    figure(1); clf
    loglog(h_step, errs_F,  'ro', 'MarkerFaceColor','r', 'MarkerSize',4); hold on
    loglog(h_step, errs_B,  'ro', 'MarkerFaceColor','g', 'MarkerSize',4);
    loglog(h_step, errs_exp_mid, 'bo', 'MarkerFaceColor','c', 'MarkerSize',4);
    loglog(h_step, errs_imp_mid, 'bo', 'MarkerFaceColor','m', 'MarkerSize',4);
    
    xlabel('Step size h');
    ylabel('Global truncation error');
    legend('Forward Euler', 'Backward', 'Explicit Midpoint', 'Implicit Midpoint')
    title('Global Truncation Error vs Step Size');

    % Error vs Step Size
    figure(2); clf
    loglog(nfe_F, errs_F,  'ro', 'MarkerFaceColor','r', 'MarkerSize',4); hold on
    loglog(nfe_B, errs_B,  'ro', 'MarkerFaceColor','g', 'MarkerSize',4);
    loglog(nfe_exp_mid, errs_exp_mid, 'bo', 'MarkerFaceColor','c', 'MarkerSize',4);
    loglog(nfe_imp_mid, errs_imp_mid, 'bo', 'MarkerFaceColor','m', 'MarkerSize',4);

    
    filter_params = struct();
    filter_params.max_xval = 1;

    [p_FE, k_FE]   = loglog_fit(h_step, errs_FE, filter_params);
    [p_Mid, k_Mid] = loglog_fit(h_step, errs_Mid, filter_params);

    loglog(h_step, k_FE*h_step.^p_FE, 'r-', 'LineWidth', 1.5);
    loglog(h_step, k_Mid*h_step.^p_Mid, 'b-', 'LineWidth', 1.5);

    legend(sprintf('Forward Euler (p = %.2f)', p_FE), sprintf('Midpoint (p = %.2f)', p_Mid),'Location','northwest');

    xlabel('function call rate');
    ylabel('global truncation error');
    legend('Forward Euler', 'Backward', 'Explicit Midpoint', 'Implicit Midpoint')
    title('Global Truncation Error vs Function Call Rate');

end