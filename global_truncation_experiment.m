
% GLOBAL: the error across the entire integration interval
% log-log fit to estimate order p


function global_truncation_experiment()
    t0 = 0;
    tf = 1;
    t_interval = [t0, tf];
    sol = @solution01;
    rate = @rate_func01;

    % Single step function methods
    forEu = @forward_euler_fixed_step_integration
    expMid = @explicit_midpoint_fixed_step_integration

    h_step = logspace(-5, 1, 100) % 1e-5 to 1e1 100 points
    figure; hold on;


        errs = zeros(size(h_step));
        for i = 1:length(h_step)
            h = h_step(i);
            XA = sol(t_interval); % finding exact X at the tref X spot
            [XB, ~] = forEu(rate, t_interval, XA, h);  % one step from exact state
            X_true = sol(t_interval + h);
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
            [XB, ~] = expMid(rate, t_interval, XA, h);  % one step from exact state
            X_true = sol(t_interval + h);
            errs(i) = norm(XB - X_true);        % norm handles vector/scalar
           
        end
        [p,k] = loglog_fit(h_step(:), errs(:));
        loglog(h_step, errs, 'k-');
        disp(p);
        disp(k);
        xlabel('h (step size)');
        ylabel('Error');
        legend('Forward Euler', 'Explicit Midpoint');

        % is y = k*x^p the same as e = O*h^p? 
        % [p,k] = loglog_fit(h_step, errs, 1)

end