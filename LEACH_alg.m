function [y, norm_Energy, CLUSTERHS,CH,GlobalMins,GlobalMin,ct,bs,out]=LEACH_alg(n,p,ETX,ERX,Efs,Emp,EDA,rmax, no_sol,dim_sol,iteration_count, S_in,objfun)
global S alg 

S=S_in;
%Computation of do
do=sqrt(Efs/Emp);

%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for rm=0:rmax
    rm
    r = rm;
    v = rm+1;

    %Operation for epoch
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
        end
    end

    hold off;

    %Number of dead nodes
    dead=0;
    %Number of dead Advanced Nodes
    dead_a=0;
    %Number of dead Normal Nodes
    dead_n=0;

    %counter for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS=0;
    packets_TO_CH=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads
    %per round
    PACKETS_TO_CH(v)=0;
    PACKETS_TO_BS(v)=0;

    figure(1);

    for i=1:1:n

        %checking if there is a dead node
        if (S(i).E<=0)
            plot(S(i).xd,S(i).yd,'red .', 'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r');
            axis off
            title(alg);
            dead=dead+1;
            if(S(i).ENERGY==1)
                dead_a=dead_a+1;
            end
            if(S(i).ENERGY==0)
                dead_n=dead_n+1;
            end
            hold on;
        end
        if S(i).E>0
            S(i).type='N';
            if (S(i).ENERGY==0)
                plot(S(i).xd,S(i).yd,'o', 'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b');
                axis off
                title(alg);
            end
            if (S(i).ENERGY==1)
                plot(S(i).xd,S(i).yd,'h', 'MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g');
                axis off
                title(alg);
            end
            hold on;
        end
    end
    plot(S(n+1).xd,S(n+1).yd,'d', 'MarkerSize',10,'MarkerEdgeColor','m','MarkerFaceColor','m');
    axis off
    title(alg);


    STATISTICS(v).DEAD=dead;
    DEAD(v)=dead;
    DEAD_N(v)=dead_n;
    DEAD_A(v)=dead_a;

    %When the first node dies
    if (dead==1)
        if(flag_first_dead==0)
            first_dead=r;
            flag_first_dead=1;
        end
    end

    countCHs=0;
    cluster=1;

    Xmax = n.*ones(no_sol,dim_sol);
    Xmin = ones(no_sol,dim_sol);
    initsol = unifrnd(Xmin,Xmax);

    %Doing process for cluster heads
    [GlobalMin(v),GlobalMins{v},GlobalParams,ct(v)]=feval(alg,initsol,objfun,Xmin,Xmax,iteration_count);
    GlobalParams = round(GlobalParams);
    GlobalParams = check_obj(n,GlobalParams);
    bs{v} = GlobalParams;
    [ObjVal,CH{v},out{v}]=feval(objfun,bs{v});

    loop=1;
    while(loop<size(GlobalParams, 2)+1)
        GlobalParams = bound_check(GlobalParams,1,n);
        i=GlobalParams(1, loop);
        loop=loop+1;
        if(S(i).E>0)

            if ( (S(i).G)<=0)
                countCHs=countCHs+1;
                packets_TO_BS=packets_TO_BS+1;
                PACKETS_TO_BS(v)=packets_TO_BS;

                S(i).type='C';
                S(i).G=round(1/p)-1;
                C(cluster).xd=S(i).xd;
                C(cluster).yd=S(i).yd;
                plot(S(i).xd,S(i).yd,'h', 'MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g');
                axis off
                %             title('alg');

                distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                C(cluster).distance=distance;
                C(cluster).id=i;
                X(cluster)=S(i).xd;
                Y(cluster)=S(i).yd;
                cluster=cluster+1;

                %Calculation of Energy dissipated
                distance;
                if (distance>do)
                    S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                end
                if (distance<=do)
                    S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                end

            end
             end
    end

    STATISTICS(v).CLUSTERHEADS=cluster-1;
    CLUSTERHS(v)=cluster-1;

    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
        if ( S(i).type=='N' && S(i).E>0 )
            if(cluster-1>=1)
                min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                min_dis_cluster=1;
                for c=1:1:cluster-1
                    temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end

                %Energy dissipated by associated Cluster Head
                min_dis;
                if (min_dis>do)
                    S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                %Energy dissipated
                if(min_dis>0)
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                    PACKETS_TO_CH(v)=n-dead-cluster+1;
                end

                S(i).min_dis=min_dis;
                S(i).min_dis_cluster=min_dis_cluster;

            end
        end
    end

    hold on;

    countCHs;
    rcountCHs=rcountCHs+countCHs;

    %Code for Voronoi Cells
    %Unfortynately if there is a small
    %number of cells, Matlab's voronoi
    %procedure has some problems

    %[vx,vy]=voronoi(X,Y);
    %plot(X,Y,'r*',vx,vy,'b-');
    % hold on;
    % voronoi(X,Y);
    % axis([0 xm 0 ym]);
    final_S=S;
    y_S=struct2cell(final_S);
    E_stru=reshape(y_S(5,1,:), [size(y_S, 3), 1]);
    E_mat=cell2mat(E_stru);
    norm_Energy(v,1)=sum(E_mat)/size(E_mat, 1);

end

x=1:1:v;
y=1:1:v;
%z=1:1:r;

for i=1:1:v
    x(i)=i;
    y(i) = n - STATISTICS(i).DEAD;
    %z(i)=CLUSTERHS(i);
end
%plot(x,y,'r',x,z,'b');
%figure, plot(x,y,'r');

end
function s = bound_check(s,LB,UB)
ns_tmp=s;
I=ns_tmp<LB;
ns_tmp(I)=LB;
J=ns_tmp>UB;
ns_tmp(J)=UB;
s=ns_tmp;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 							  %
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round					  %
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round                     %
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round                      %
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round      %
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round   %
%  first_dead: the round where the first node died                                    %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%