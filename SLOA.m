function [best_position, best_score, Convergence_curve, time] = SLOA(population, objective_function, lb, ub, max_iter)
% Snow Leopard Optimization Algorithm (SLOA)
% objective_function: function handle of the objective function
% num_agents: number of snow leopards
% max_iter: maximum number of iterations
% lb: lower bound of variables (vector)
% ub: upper bound of variables (vector)
% dim: number of variables (dimension of the problem)
[num_agents,dim] = size(population);
fitness = zeros(num_agents, 1);
for i = 1:num_agents
    fitness(i) = feval(objective_function, population(i, :));
end

% Initialize the best solution
[best_score, best_idx] = min(fitness);
best_position = population(best_idx, :);
Convergence_curve=zeros(1, max_iter);
tic;
% Main loop of the SLOA
for iter = 1:max_iter
    for i = 1:num_agents
        % Phase 1: Travel Routes and Movement
        r = rand();
        k = randi(num_agents);
        I = round(1 + r);
        new_position = population(i, :) + r * (population(k, :) - I * population(i, :)) .* sign(fitness(i) - fitness(k));
        new_position = bound_check(new_position, lb, ub);
        new_fitness = feval(objective_function, new_position);
        if new_fitness < fitness(i)
            population(i, :) = new_position;
            fitness(i) = new_fitness;
        end
        
        % Phase 2: Hunting
        prey_position = best_position;
        P = 0.375; % Walking proportion
        Q = 0.625; % Running proportion
        for d = 1:dim
            r1 = rand();
            if r1 < P
                new_position(d) = population(i, d) + r * (prey_position(d) - population(i, d));
            else
                new_position(d) = prey_position(d) - r * (prey_position(d) - population(i, d));
            end
        end
        new_position = bound_check(new_position, lb, ub);
        new_fitness = feval(objective_function, new_position);
        if new_fitness < fitness(i)
            population(i, :) = new_position;
            fitness(i) = new_fitness;
        end
        
        % Phase 3: Reproduction
        mate_idx = randi(num_agents);
        offspring = 0.5 * (population(i, :) + population(mate_idx, :));
        offspring = bound_check(offspring, lb, ub);
        offspring_fitness = feval(objective_function, offspring);
        if offspring_fitness < fitness(i)
            population(i, :) = offspring;
            fitness(i) = offspring_fitness;
        end
        
        % Phase 4: Mortality
        if fitness(i) > mean(fitness)
            population(i, :) = lb + rand(1, dim) .* (ub - lb);
            fitness(i) = feval(objective_function, population(i, :));
        end
    end
    
    % Update the best solution
    [current_best_score, best_idx] = min(fitness);
    if current_best_score < best_score
        best_score = current_best_score;
        best_position = population(best_idx, :);
    end
    
    % Display iteration information
    disp(['Iteration ' num2str(iter) ': Best Score = ' num2str(best_score)]);
end
iter=iter+1;
Convergence_curve(iter)=best_score;

time = toc;
end

function position = bound_check(position, lb, ub)
% Ensure the position is within the bounds
position = max(position, lb);
position = min(position, ub);
end
