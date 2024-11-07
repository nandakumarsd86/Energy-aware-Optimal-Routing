function [best_fitness, best_solution, Convergence_curve, time] = PROPOSED(X, sphere_function, lb, ub, max_iterations)
% Dark Forest Algorithm (DFA)
[n, dim] = size(X);
% Initialize parameters
b = 1; % Logarithmic spiral shape constant
delta = 0.5; % Mapping factor
beta = 0.5; % Search radius
iter = 0;

% Compute initial fitness
fitness = zeros(n, 1);
for i = 1:n
    fitness(i) = feval(sphere_function, X(i, :));
end

[sorted_fitness, sorted_indices] = sort(fitness);
X = X(sorted_indices, :);
Convergence_curve = zeros(1, max_iterations);
tic;
stepsize=zeros(n,dim);
% Main loop
while iter < max_iterations - 10
    % Classification of civilizations
    highest = X(1, :);
    advanced = X(2:ceil(0.3*n), :);
    normal = X(ceil(0.3*n)+1:ceil(0.8*n), :);
    low = X(ceil(0.8*n)+1:end, :);
    
    % Update positions
    % Highest civilization
    for i = 1:size(highest, 1)
        if rand < 0.5
            highest(i, :) = highest(i, :) + randn(1, dim);
        else
            ref_advanced = advanced(randi(size(advanced, 1)), :);
            highest(i, :) = highest(i, :) + rand(1, dim) .* (ref_advanced - highest(i, :));
        end
        highest(i, :) = min(max(highest(i, :), lb), ub);
    end
    %
    r1 = min(fitness) / (max(fitness) ^ 2);
    r2 = min(fitness) / (max(fitness) + 2);
    r3 = r1 / (r2 ^ 3);
    
    
    Elite=repmat(X,n,1);  %(Eq. 3)
    RB=randn(n,dim);          %Brownian random number vector
    % Advanced civilizations
    for i = 1:size(advanced, 1)
        if r3 < r2
            ref_normal = normal(randi(size(normal, 1)), :);
            advanced(i, :) = ref_normal;
        else
            % update using GOA
            
            for j=1:size(gazelle,2)
                stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*advanced(i,j));
                advanced(i,j)=advanced(i,j)+s*R*stepsize(i,j);
            end
        end
        advanced(i, :) = min(max(advanced(i, :), lb), ub);
    end
    
    % Normal civilizations
    for i = 1:size(normal, 1)
        ref_civ = X(randi(size(X, 1)), :);
        D = norm(ref_civ - normal(i, :));
        normal(i, :) = normal(i, :) + D * exp(b * rand * (2 * pi)) .* cos(2 * pi * rand(1, dim));
        normal(i, :) = min(max(normal(i, :), lb), ub);
    end
    
    % Low civilizations
    for i = 1:size(low, 1)
        if rand < 0.5
            low(i, :) = lb + (ub - lb) .* rand(1, dim);
        else
            Xw = mean([highest; advanced], 1);
            low(i, :) = Xw + delta * (Xw - low(i, :)) .* (1 - iter/max_iterations)^2;
        end
        low(i, :) = min(max(low(i, :), lb), ub);
    end
    
    % Combine all civilizations
    X = [highest; advanced; normal; low];
    
    % Compute fitness
    for i = 1:n
        fitness(i) = feval(sphere_function, X(i, :));
    end
    
    [sorted_fitness, sorted_indices] = sort(fitness);
    X = X(sorted_indices, :);
    
    % Update best solution
    best_solution = X(1, :);
    best_fitness = sorted_fitness(1);
    Convergence_curve(iter+1) = best_fitness;
    fprintf('Iteration %d: Best Fitness = %f\n', iter+1, best_fitness);
    
    iter = iter + 1;
end

% Refine search for the highest civilization
for k = 1:10
    for i = 1:dim
        beta_i = beta;
        X_new = best_solution;
        X_new(i) = best_solution(i) + beta_i;
        if feval(sphere_function, X_new) < best_fitness
            best_solution = X_new;
            best_fitness = feval(sphere_function, X_new);
        else
            beta_i = -beta_i;
            X_new(i) = best_solution(i) + beta_i;
            if feval(sphere_function, X_new) < best_fitness
                best_solution = X_new;
                best_fitness = feval(sphere_function, X_new);
            else
                beta_i = beta_i / 2;
            end
        end
    end
end

time = toc;
end
