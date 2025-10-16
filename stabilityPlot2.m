function stabilityPlot2()

    tspan = [0, 10];
    X0 = 1;
    h = 0.45;
    
    forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h)
    explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h)

    fixed_step_integration(@rate_func01, @backward_euler_step,tspan,X0,h)
    fixed_step_integration(@rate_func01, @implicit_midpoint_step,tspan,X0,h)

    
    t_exact = linspace(tspan(1), tspan(2), 1000); % Create a vector of 1000 evenly spaced points in the interval [tspan(1),tspan(12)].

    X_exact = solution01(t_exact);

    figure(1);
    plot(t_exact, X_exact, 'k-', 'LineWidth', 1.5); hold on;
    
    [t_forEu, X_forEu] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h);
    plot(t_forEu, X_forEu,'g--', 'LineWidth', 1.5); hold on;

    legend('Exact', 'h-ref = 0.45');
    xlabel('Time t'); ylabel('x(t)');
    title('Forward Euler');


    figure(2);
    plot(t_exact, X_exact, 'k-', 'LineWidth', 1.5); hold on;
 
    [t_expMid, X_expMid] =  explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0, h);

    plot(t_expMid, X_expMid, 'r--', 'LineWidth', 1.5); hold on;
  
    legend('Exact', 'h-ref = 0.45');
    xlabel('Time t'); ylabel('x(t)');
    title('Explicit Midpoint');

    
    figure(3);
    plot(t_exact, X_exact,'k-', 'LineWidth', 1.5); hold on;
     
    [t_impEu, X_impEu, h_avg, num_evals] = fixed_step_integration(@rate_func01,@backward_euler_step, tspan,X0,h);
    plot(t_impEu, X_impEu,'b--', 'LineWidth', 1.5); hold on;

    legend('Exact', 'h-ref = 0.45');
    xlabel('Time t'); ylabel('x(t)');
    title('Backward Euler');

    
    figure(4);
    plot(t_exact, X_exact, 'k-', 'LineWidth', 1.5); hold on;

    [t_impMid, X_impMid, h_avg, num_evals] = fixed_step_integration(@rate_func01,@backward_euler_step, tspan,X0,h);
    plot(t_impMid, X_impMid, 'y--', 'LineWidth', 1.5); hold on;

    legend('Exact', 'h-ref = 0.45');
    xlabel('Time t'); ylabel('x(t)');
    title('Implicit Midpoint');


    % To verify that your implementation is working, 
    % use it to solve the IVP described by equation 5. 
    % Plot the closed-form solution for X(t) 
    % and a few numerical approximations 
    % (with different time steps) on the same axes.
end