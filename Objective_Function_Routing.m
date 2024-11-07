function [ObjVal, cluster_center, out] = Objective_Function_Routing(soln)
global S n Source Dest ResidualEnergy RouteQuality RouteEnergy Congestion SINR

for col = 1:size(soln, 1)
    CH = round(soln(col, :));
    CH = check_obj(n, CH);
    CH = bound_check(CH, 1, n);
    Path = unique(CH);
    short_path = [Source, Path, Dest];  % the path with source and dest
    
    xx = S;
    y2 = struct2cell(xx);
    xd = reshape(y2(1, 1, :), [size(y2, 3), 1]);
    yd = reshape(y2(2, 1, :), [size(y2, 3), 1]);
    G = reshape(y2(3, 1, :), [size(y2, 3), 1]);
    type = reshape(y2(4, 1, :), [size(y2, 3), 1]);
    E = reshape(y2(5, 1, :), [size(y2, 3), 1]);
    Energy = reshape(y2(6, 1, :), [size(y2, 3), 1]);
    cl = reshape(y2(7, 1, :), [size(y2, 3), 1]);
    
    loca_data = [cell2mat(xd) cell2mat(yd)];
    cluster_center = loca_data(CH, :);
    
    %% Distance Constraints
    distance_Ch_node = zeros(size(y2, 3) - 1, size(CH, 2));
    for i = 1:size(loca_data, 1) - 1
        for j = 1:size(CH, 2)
            distance_Ch_node(i, j) = dist(cluster_center(j, :), loca_data(i, :)');
        end
    end
    
    for j = 1:size(CH, 2)
        distance_CH_base(j, 1) = dist(cluster_center(j, :), loca_data(end, :)');
    end
    
    a = cluster_center(:, 1);
    b = cluster_center(:, 2);
    for l = 1:size(a, 1)
        for m = 1:size(b, 1)
            distn(m) = dist(a(l), b(m));
        end
        Idist(l) = sum(distn);
    end
    Idist = sum(Idist) / 10;
    
    %% Finding fitness-distance
    [min_value, index_val] = min(distance_Ch_node, [], 2);
    count = zeros(size(CH, 2), 1);
    
    for i = 1:size(loca_data, 1) - 1
        for j = 1:size(CH, 2)
            if j == index_val(i)
                f_dist_b = distance_CH_base(j, 1) + distance_Ch_node(i, j);
                count(j, 1) = count(j, 1) + 1;
            end
        end
    end
    
    f_dist = 1 / f_dist_b;
    
    %% Finding Energy
    E_mat = cell2mat(E);
    CH_E = E_mat(CH, :);
    [E_CH_min_value, E_CH_index_val] = min(CH_E, [], 1);
    
    [E_node_min_value, E_node_index_val] = min(E_mat, [], 1);
    phi = 20.72;
    
    tou_value = -phi * (E_CH_min_value / (abs(E_node_min_value - E_CH_min_value) + 1e-10));
    f_energy_b = 1 - exp(tou_value);
    
    f_energy = abs(f_energy_b / exp(sum(count)));
    %% Residual Energy
    ResidualEnergyVal = sum(ResidualEnergy(short_path));
    
    %% Route Quality
    RouteQualityVal = sum(RouteQuality(short_path));
    
    %% Route Energy
    RouteEnergyVal = sum(RouteEnergy(short_path));
    
    %% Congestion
    CongestionVal = sum(Congestion(short_path));
    
    %% SINR
    SINRVal = sum(SINR(short_path));
    
    %% Applying normalization before unification
    alpha = 0.2;
    beta = 0.2;
    gamma = 0.2;
    delta = 0.2;
    omega = 0.2;
    F1 = alpha * (1 / ResidualEnergyVal) + (1 - alpha) * (1 / ResidualEnergyVal);
    F2 = beta * F1 + (1 - beta) * (1 / RouteQualityVal);
    F3 = gamma * F2 + (1 - gamma) * (RouteEnergyVal);
    F4 = delta * F3 + (1 - delta) * (CongestionVal);
    F5 = omega * F4 + (1 - omega) * (SINRVal);
    ObjVal(col) = 1 / F5;
    
    out = [ResidualEnergyVal, RouteQualityVal, RouteEnergyVal, CongestionVal, SINRVal];
end
end

function s = bound_check(s, LB, UB)
ns_tmp = s;
ns_tmp(isnan(ns_tmp)) = 1;
I = ns_tmp < LB;
ns_tmp(I) = LB;
J = ns_tmp > UB;
ns_tmp(J) = UB;
s = ns_tmp;
end

function d = dist(x1, x2)
d = sqrt(sum((x1 - x2).^2));
end
