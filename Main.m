function [] = Main_17_05_2023()
%% Main Function
clc;
clear all;
close all;
warning off

global n Source Dest Total_Packets Packet_loss  ETX Efs Emp ERX S alg

num_of_node = [50, 100, 150, 200];

%% Experiments
an = 0;
if an == 1
    %Field Dimensions - x and y maximum (in meters)
    xm=100;
    ym=100;
    
    %maximum number of rounds
    rmax = 2000;
    
    %Optimization paramateres
    no_sol=10;
    dim_sol=30;
    iteration_count=10;
    
    %x and y Coordinates of the Sink
    sink.x=0.5*xm;
    sink.y=0.5*ym;
    
    for net = 1:length(num_of_node)
        %Number of Nodes in the field
        n = num_of_node(net);
        
        %Optimal Election Probability of a node
        %to become cluster head
        p=0.1;
        
        Total_Packets = 1000;
        Packet_loss = randi([0, 5], 1, n);
        
        %Energy Model (all values in Joules)
        %Initial Energy
        Eo=0.3;
        
        %Eelec=Etx=Erx
        ETX=50*0.000000001;
        ERX=50*0.000000001;
        %Transmit Amplifier types
        Efs=10*0.000000000001;
        Emp=0.0013*0.000000000001;
        %Data Aggregation Energy
        EDA=5*0.000000001;
        
        %Values for Hetereogeneity
        %Percentage of nodes than are advanced
        m=0.1;
        %alpha
        a=1;
        
        %% Creation of the random Sensor Network
        figure(1);
        for i=1:1:n
            S_in(i).xd=rand(1,1)*xm;
            XR(i)=S_in(i).xd;
            S_in(i).yd=rand(1,1)*ym;
            YR(i)=S_in(i).yd;
            S_in(i).G=0;
            %initially there are no cluster heads only nodes
            S_in(i).type='N';
            
            temp_rnd0=i;
            %Random Election of Normal Nodes
            if (temp_rnd0>=m*n+1)
                S_in(i).E=Eo;
                S_in(i).ENERGY=0;
                plot(S_in(i).xd,S_in(i).yd,'o', 'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b');
                axis off
                hold on;
            end
            %Random Election of Advanced Nodes
            if (temp_rnd0<m*n+1)
                S_in(i).E=Eo*(1+a);
                S_in(i).ENERGY=1;
                plot(S_in(i).xd,S_in(i).yd,'h', 'MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g');
                axis off
                hold on;
                
            end
        end
        S_in(n+1).xd=sink.x;
        S_in(n+1).yd=sink.y;
        plot(S_in(n+1).xd,S_in(n+1).yd,'d', 'MarkerSize',10,'MarkerEdgeColor','m','MarkerFaceColor','m');
        axis off
        
        %% Cluster formation
        xx=S_in;%load('S.mat');
        y2=struct2cell(xx);
        xd=reshape(y2(1,1,:), [size(y2, 3), 1]);
        yd=reshape(y2(2,1,:), [size(y2, 3), 1]);
        loca_data=[cell2mat(xd) cell2mat(yd)];
        [idx,C] = kmeans(loca_data,10);
        for j = 1:length(C)
            A = repmat(C(j,:),length(loca_data),1);
            [minValue,CH_indd(j)] = min(sum(abs(loca_data-A)'));
            CHs(j,:) = loca_data(CH_indd(j),:);  % Cluster Heads
            S_in(CH_indd(j)).E=0.5;  % Imitial Energy of CH is high than ordinary node
        end
        S_in(n+1).xd=sink.x;
        S = S_in;
        %% Analysis
        %Optimization for  Routing
        an = 0;
        if an ==1
            Algms = {'MOA', 'SLOA', 'DFA', 'GOA', 'PROPOSED'}; 
            for i = 1:length(Algms)
                alg = Algms{i};
                MainNodes = randperm(n, 2);
                Source = MainNodes(1);
                Dest = MainNodes(2);
                [y_An, norm_Energy, CLUSTERHS,CH,GM,F,ct,bs, out]=LEACH_alg(n,p,ETX,ERX,Efs,Emp,EDA,rmax, no_sol,dim_sol,iteration_count,S_in,'Objective_Function_Routing');
                Result_ct(net,i).y_An = y_An;
                Result_ct(net,i).norm_Energy = norm_Energy;
                Result_ct(net,i).CLUSTERHS = CLUSTERHS;
                Result_ct(net,i).CH = CH;
                Result_ct(net,i).GM = GM;
                Result_ct(net,i).F = F;
                Result_ct(net,i).tm = ct;
                Result_ct(net,i).sol = bs;
                Result_ct(net,i).S = S;
                Result_ct(net,i).out = out;
                save Result_Ch Result_ct
            end
        end
    end
end
Plot_Results()

end

