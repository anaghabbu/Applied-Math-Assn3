%This function computes the value of X at the next time step
%using the Backward Euler approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)

    G = @(X_in) XA +  h*rate_func_in(t + h, X_in)  -  X_in;
    
    solver_params = struct();
    solver_params.dxtol = 1e-6;
    solver_params.ftol = 1e-6;
    solver_params.max_iter = 50;
    solver_params.dxmax = 1e3;
    solver_params.numerical_diff = true;

    [XB, num_evals]  = multi_newton_solver2(G,XA,solver_params);    

    % G_xnplus1 = XA + h*rate_func_in(t + h, )

  
    % G(Xn+1) = Xn + hf (tn + h, Xn+1) − Xn+1 = 0
end