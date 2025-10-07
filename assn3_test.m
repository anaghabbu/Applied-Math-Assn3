function assn3_test()

    tspan = [3, 12];
    X0 = 1;
    h_ref = 2;

    forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref)

    explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref)

end