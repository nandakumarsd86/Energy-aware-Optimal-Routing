function [best_fitness, best_solution,Convergence_curve, time]=MOA(mothers,sphere_function,lb,ub,max_iterations)
% Mother Optimization Algorithm (MOA)

[num_mothers,num_variables] = size(mothers);
% Parameters
num_children = 5; % Number of children per mother
mutation_rate = 0.1; % Mutation rate
% Initialize population
children = zeros(num_mothers * num_children, num_variables); % Initialize children

Convergence_curve=zeros(1, max_iterations);
tic;
% Main loop
for iteration = 1:max_iterations
    % Generate children
    for i = 1:num_mothers
        for j = 1:num_children
            child_index = (i - 1) * num_children + j;
            mother = mothers(i, :);
            child = mother + mutation_rate * randn(1, num_variables); % Mutation
            % Ensure children stay within bounds
            child = min(max(child, lb), ub);
            children(child_index, :) = child(child_index);
        end
    end
    
    % Evaluate objective function for children
    child_fitness = zeros(size(children, 1), 1);
    for i = 1:size(children, 1)
        child_fitness(i) = feval(sphere_function, children(i, :));
    end
    
    % Combine mothers and children
    combined_population = [mothers; children];
    
    % Sort combined population based on fitness
    [~, sorted_indices] = sort(child_fitness);
    combined_population = combined_population(sorted_indices, :);
    
    % Select best mothers
    mothers = combined_population(1:num_mothers, :);
    
    % Display best solution
    best_solution = mothers(1, :);
    best_fitness = feval(sphere_function, best_solution);
    fprintf('Iteration %d: Best Fitness = %f\n', iteration, best_fitness);
end

iteration=iteration+1;
Convergence_curve(iteration)=best_fitness;

time = toc;
end