function explicit_v_implicit_trunc()

    tref = 0.492;
    sol = @solution01;
    rate = @rate_func01;

    % Single step function methods
    forward_Eu = @forward_euler_step
    backward_Eu = @backward_euler_step

    exp_mid = @explicit_midpoint_step
    imp_mid = @implicit_midpoint_step

    h_step = logspace(-5, 1, 100); % 1e-5 to 1e1 100 points
    
    %analytical_diff  = zeros(size(h_step));
    errs_forward = zeros(size(h_step));
    errs_backward = zeros(size(h_step));
    errs_exp_mid = zeros(size(h_step));
    errs_imp_mid = zeros(size(h_step));

    XA = sol(tref); % finding exact X at the tref X spot

    for r = 1:length(h_step)
        h = h_step(r);

        
        [XB_F, ~] = forward_Eu(rate, tref, XA, h);  % one step from exact state
        [XB_B, ~] = backward_Eu(rate, tref, XA, h);  % one step from exact state
        [XB_x, ~] = exp_mid(rate, tref, XA, h);  % one step from exact state
        [XB_i, ~] = imp_mid(rate, tref, XA, h);  % one step from exact state
        X_true = sol(tref + h);

        %analytical_diff(r) = norm(X_true - XA);
        errs_forward(r) = norm(XB_F - X_true);
        errs_backward(r) = norm(XB_B - X_true);
        errs_exp_mid(r) = norm(XB_x - X_true);
        errs_imp_mid(r) = norm(XB_i - X_true);


        % errs(i) = norm(XB - X_true);        % norm handles vector/scalar

    end
    % [p,k] = loglog_fit(h, errs(i), 1);
    % plot(errs(i), h);
    
    figure(1); 
    %loglog(h_step,analytical_diff,'ko','MarkerFaceColor','k','MarkerSize',1);
    %hold on
    loglog(h_step,errs_forward,'r-','MarkerFaceColor','r','MarkerSize',1); 
    hold on
    loglog(h_step,errs_backward,'b-','MarkerFaceColor','b','MarkerSize',1);

    %filter_params = struct();
    %filter_params.max_xval = 1;

    %figure(2);
    loglog(h_step,errs_exp_mid,'g-','MarkerFaceColor','r','MarkerSize',1); 
    %hold on
    loglog(h_step,errs_imp_mid,'m-','MarkerFaceColor','b','MarkerSize',1);

    filter_params = struct();
    filter_params.max_xval = 1;

    legend('Forward Euler','Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint')
    title('Local Truncation Error for Explicit and Implicit Methods')

    %[p1,k1] = loglog_fit(h_step, analytical_diff, filter_params);
    %[p2,k2] = loglog_fit(h_step, errs_forward, filter_params);
    %[p3,k3] = loglog_fit(h_step, errs_backward, filter_params);

    %loglog(h_step,k1*h_step.^p1,'g','LineWidth',1)
    %loglog(h_step,k2*h_step.^p2,'r','LineWidth',1)
    %loglog(h_step,k3*h_step.^p3,'b','LineWidth',1)


    % is y = k*x^p the same as e = O*h^p? 
    % [p,k] = loglog_fit(h_step, errs, 1)
 

end