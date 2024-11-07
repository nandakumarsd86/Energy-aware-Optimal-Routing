function [] = Plot_Results()

Plot_Cluster_results()
Convergence_Graph()

end

%% Results for Routing
function [] = Plot_Cluster_results()
load Results_ch

Stats = {'BEST', 'WORST', 'MEAN', 'MEDIAN', 'STD'};
Terms = {'Energy Consumption (j)' 'Delay' 'Residual Energy' 'Throughput (%)' 'Route Energy' 'Route quality' 'congestion' 'SINR' 'Packet Lost' 'Packets Obtained' 'packet delivery ratio (%)'};
num_of_nodes = [50, 100, 150, 200];

for n = 1:4
    for i = 1:size(Results_ch, n)
        for j = 1:length(Results_ch(n, i).out)
            Outs{i}(j, :) = Results_ch(n, i).out{j};
        end
        Outs{i}(isinf(Outs{i})|isnan(Outs{i})) = 0;
    end
end

Positive = [3, 4, 5,6,11];

for i = 1:length(Terms)
    for j = 1:length(Outs)
        if length(find(ismember(i, Positive))) >= 1
            Statistics(j, 1) = max(Outs{j}(:, i));
            Statistics(j, 2) = min(Outs{j}(:, i));
            Statistics(j, 3) = mean(Outs{j}(:, i));
            Statistics(j, 4) = median(Outs{j}(:, i));
            Statistics(j, 5) = std(Outs{j}(:, i));
        else
            Statistics(j, 1) = min(Outs{j}(:, i));
            Statistics(j, 2) = max(Outs{j}(:, i));
            Statistics(j, 3) = mean(Outs{j}(:, i));
            Statistics(j, 4) = median(Outs{j}(:, i));
            Statistics(j, 5) = std(Outs{j}(:, i));
        end
    end
    disp(strcat("-------------------- Number of node - ", num2str(num_of_nodes(n)), " - ", char(Terms{i}), " Statistical Report --------------------"))
    T = table(char(Stats), Statistics(1, :)', Statistics(2, :)', Statistics(3, :)', Statistics(4, :)', Statistics(5, :)');
    T.Properties.VariableNames = {'Statistics', 'MOA', 'SLOA', 'DFA', 'GOA', 'Proposed'};
    disp(T)
end



X = 50:50:200;
for i = 1:length(Terms)
    for j = 1:size(Results_ch, 2)
        Plt_Vaues(j, :) = sort(Outs{j}(X, i));
        
    end
    for j = 1:size(Plt_Vaues, 2)
        if length(find(ismember(i, Positive))) >= 1
            [a, b] = max(Plt_Vaues(:, j));
        else
            [a, b] = min(Plt_Vaues(:, j));
        end
        x = Plt_Vaues(b, j);
        y = Plt_Vaues(end, j);
        Plt_Vaues(end, j) = x;
        Plt_Vaues(b, j) = y;
    end
    figure,
    bar(X, Plt_Vaues')
    newColors = [0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.75, 0.75, 0; 1, 0, 0; 0.75, 0, 0.75];
    colororder(newColors)
    set(gca, 'FontSize', 14);
    xlabel('Number of nodes', 'FontSize', 14);
    ylabel(char(Terms{i}), 'FontSize', 14);
    h = legend('MOA', 'SLOA', 'DFA', 'GOA', 'Proposed');
    set(h, 'fontsize', 12, 'Location', 'NorthEastOutside')
    print('-dtiff', '-r300', ['.\Results\', char(Terms{i}), '-', num2str(n)])
end
end


function [] = Convergence_Graph()
Node  = [50,100,150,200];
Dataset = [];
load fitness1;
for n = 1 : size(Fit,2)
    for j = 1 : size(Fit{1, 1},1) % For all algorithms
        val(j,:) = stats(Fit{1, n}(j,:));
    end
    disp('Statistical Analysis :')
    fprintf('Number of Node  : %d\n ', Node (n));
    ln = {'BEST','WORST','MEAN','MEDIAN','STANDARD DEVIATION'};
    T = table(val(1, :)', val(2, :)', val(3, :)',val(4, :)', val(5, :)','Rownames',ln);
    T.Properties.VariableNames = {'BMO', 'OOA', 'DFA', 'GOA', 'Proposed'};
    disp(T)

    figure,
    plot(Fit{1, n}(1,:),'m', 'markersize', 5, 'LineWidth', 4)
    hold on;
    plot(Fit{1, n}(2,:),'g', 'markersize', 5, 'LineWidth', 4)
    plot(Fit{1, n}(3,:),'r', 'markersize', 5, 'LineWidth', 4)
    plot(Fit{1, n}(4,:),'b', 'markersize', 5, 'LineWidth', 4)
    plot(Fit{1, n}(5,:),'k', 'markersize', 5, 'LineWidth', 4)
    set(gca,'fontsize',20);
    xlabel('No. of Iterations','fontsize',16);
    xticks([20 40 60 80 100])
    ylabel('Cost Function','fontsize',16);
    h = legend('BMO', 'OOA', 'DFA', 'GOA', 'Proposed');
    set(h,'fontsize',12,'Location','Best')
    print('-dtiff','-r300',['.\Results\', 'Convergence-',num2str(n)])
end
end


function[a] = stats(val)
a(1) = min(val);
a(2) = max(val);
a(3) = mean(val);
a(4) = median(val);
a(5) = std(val);
end


