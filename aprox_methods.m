% For Galerkin Method (l = 1)
aprox_methods(1, 1)

% For Least Squares Method (l = 2)
aprox_methods(1, 2)

% For Collocation Method (l = 3)
aprox_methods(1, 3)


function u_at_x_05 = aprox_methods(N, l)
    % Function that implements Galerkin, Least Squares, and Collocation methods
    % N: number of terms in the series
    % l: type of method (1 = Galerkin, 2 = Least Squares, 3 = Collocation)
    
    % Defining symbolic variables
    syms x a;

    % Initializing the approximator function u(x)
    u = 0; 

    % Loop to add terms in the series (summing up to the N-th term)
    for k = 1:N
        u = u + a^k * x^k * (1 - x); % Adds the term a_k * x^k * (1 - x)
    end

    % Defining the equation representing the residue f
    f = diff(u, x, 2) - u + x; % Second derivative of u minus u plus x

    % Weight function w, which depends on the chosen method
    switch l
        case 1  % Galerkin Method
            w = diff(u, a);  % Derivative of u with respect to a
        case 2  % Least Squares Method
            w = diff(f, a);  % Derivative of f with respect to a
        case 3  % Collocation Method
            % Defining the collocation point, for example x = 0.5
            x_i = 0.5;
            
            % Evaluating the residue at x_i
            residuo = subs(f, x, x_i);  % Substituting x_i in the residue equation
            
            % The residue must be forced to zero at x_i
            Aprox = residuo;  % Residue at x_i
    end

    % For Galerkin and Least Squares methods, we use the integral.
    if l ~= 3
        % Multiplying the residue by the weight function
        Aprox = f * w;

        % Defining the function interval from 0 < x < 1
        inter1 = 0;
        inter2 = 1;

        % Defining the values of the approximator function (integrating the residue times the weight function)
        L = int(Aprox, x, inter1, inter2);

        % Displaying the integral value (solving for L = 0)
        a_n = solve(L == 0, a);  % Solving the integral to find the solution for a
    else
        % In the Collocation method, integration is not needed
        a_n = solve(Aprox == 0, a);  % Solving for a_n
    end

    % Displaying the value of a_n (solution for a)
    disp('Value of a_n (solution for a):');
    disp(a_n);  

    % Substituting the value of a_n into the function u(x)
    u_sub = subs(u, a, a_n);  % Substituting a_n into u

    % Evaluating the approximator function at x = 0.5
    u_at_x_05 = double(subs(u_sub, x, 0.5));  % Evaluating the function at x = 0.5

    % Displaying the value of u at x = 0.5
    disp('Value of u at x = 0.5:');
    disp(u_at_x_05);  
end
