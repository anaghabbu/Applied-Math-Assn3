% To verify that your implementation is working, use it to solve the IVP described by equation 5. Plot the closed-form
% solution for X(t) and a few approximations (with different time steps) on the same axes.

function assn3_test()

    tspan = [0, 10];
    X0 = 1;

    h_values = [0.2, 0.1, 0.05, 0.02]; % using a few diff time steps

    forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref)
    explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref)

    fixed_step_integration(@rate_func01, @backward_euler_step,tspan,X0,h_ref)
    fixed_step_integration(@rate_func01, @implicit_midpoint_step,tspan,X0,h_ref)

    
    t_exact = linspace(tspan(1), tspan(2), 1000); % Create a vector of 1000 evenly spaced points in the interval [tspan(1),tspan(12)].

    X_exact = solution01(t_exact);

    figure(1);
    plot(t_exact, X_exact, 'k-', 'LineWidth', 1.5); hold on;
    randColor = rand(10,3);
    for i = 1:length(h_values)
        h = h_values(i);
        [t_forEu, X_forEu] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h);
        plot(t_forEu, X_forEu,'Color',randColor(i,:), 'LineWidth', 1.5); hold on;
    end
    legend('Exact', 'Timestep 0.2', 'Timestep 0.1', 'Timestep 0.05', 'Timestep 0.02');
    xlabel('Time t'); ylabel('x(t)');
    title('Forward Euler');

    figure(2)
    plot(t_exact, X_exact, 'k-', 'LineWidth', 1.5); hold on;
    randColor = rand(10,3);
    for i = 1:length(h_values)
        h = h_values(i);
        [t_expMid, X_expMid] =  explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0, h);
        figure(2);
        plot(t_expMid, X_expMid, 'Color',randColor(i,:),'LineWidth', 1.5); hold on;

    end
    legend('Exact', 'Timestep 0.2', 'Timestep 0.1', 'Timestep 0.05', 'Timestep 0.02');
    xlabel('Time t'); ylabel('x(t)');
    title('Explicit Midpoint');

    figure(3)
    plot(t_exact, X_exact,'Color',randColor(i,:), 'LineWidth', 1.5); hold on;
    randColor = rand(10,3);
    for i = 1:length(h_values)
        h = h_values(i);
        
        [t_impEu, X_impEu, h_avg, num_evals] = fixed_step_integration(@rate_func01,@backward_euler_step, tspan,X0,h);
      
        figure(3);
        plot(t_impEu, X_impEu,'Color',randColor(i,:), 'LineWidth', 1.5); hold on;


    end
    legend('Exact', 'Timestep 0.2', 'Timestep 0.1', 'Timestep 0.05', 'Timestep 0.02');
    xlabel('Time t'); ylabel('x(t)');
    title('Backward Euler');

    figure(4)
    plot(t_exact, X_exact, 'k-', 'LineWidth', 1.5); hold on;
    randColor = rand(10,3);
    for i = 1:length(h_values)
        h = h_values(i);
        
        [t_impMid, X_impMid, h_avg, num_evals] = fixed_step_integration(@rate_func01,@backward_euler_step, tspan,X0,h);
     
        figure(4)
      
        plot(t_impMid, X_impMid, 'Color',randColor(i,:), 'LineWidth', 1.5); hold on;

    end
    legend('Exact', 'Timestep 0.2', 'Timestep 0.1', 'Timestep 0.05', 'Timestep 0.02');
    xlabel('Time t'); ylabel('x(t)');
    title('Implicit Midpoint');


    % To verify that your implementation is working, 
    % use it to solve the IVP described by equation 5. 
    % Plot the closed-form solution for X(t) 
    % and a few numerical approximations 
    % (with different time steps) on the same axes.
end