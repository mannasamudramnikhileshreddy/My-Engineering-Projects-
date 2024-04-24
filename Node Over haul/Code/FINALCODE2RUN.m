
close all;
clear all;
xm=50;      %diameters of sensor network
ym=50;
sink.x=10;  %distance of base station from the network
sink.y=10;
n = 50;         %no of nodes
p=0.1;          %probibilty of a node to become cluster head
Eo=0.28;          %energy supplied to each node
ETX=50*0.000000001;     %transmiter energy per node
ERX=50*0.000000001;        %reciever energy per mode
Efs=10*0.0000000000001;     %amplification energy when d is less than d0
Emp=0.0013*0.000000000001;      %amplification energy  when d is greater than d0
%Data Aggregation Energy
EDA=5*0.000000001;
rmax=1200;           %no of rounds
do=sqrt(Efs/Emp);       %distance between cluster head and base station

for i=1:1:n
    S(i).xd=rand(1,1)*xm;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S(i).yd=rand(1,1)*ym;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S(i).G=0;                        % as the no of node that have been cluster head is zero 0
    S(i).E=Eo;
    %initially there are no cluster heads only nodes
    S(i).type='N';
end
S(n+1).xd=sink.x;   %assume that base station is also a node sp total no of nodes is n and with base station  it is n+1
S(n+1).yd=sink.y;
countCHs=0;         %the number of Stateflow objects in the current context.
cluster=1;              %first cluster is selected
flag_first_dead=0;
flag_teenth_dead=0;
flag_all_dead=0;
dead=0;
first_dead=0;
teenth_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
% counter for sleep nodes
s=0;

for r=0:1:rmax
    r
    if(mod(r, round(1/p) )==0) %remainder
        for i=1:1:n
            S(i).G=0;            % it will assign to the nodes that have not been cluster head .
        end
    end
    dead=0;
    for i=1:1:n
        
        if (S(i).E<=0)
            dead=dead+1;
            
            if (dead==1)
                if(flag_first_dead==0)
                    first_dead=r;
                    flag_first_dead=1;
                end
            end
            
            if(dead==0.1*n)
                if(flag_teenth_dead==0)
                    teenth_dead=r;
                    flag_teenth_dead=1;
                end
            end
            if(dead==n)
                if(flag_all_dead==0)
                    all_dead=r;
                    flag_all_dead=1;
                end
            end
        end
        if (S(i).E>0)
            S(i).type='N';
        end
    end
    STATISTICS.DEAD(r+1)=dead;
    STATISTICS.ALLIVE(r+1)=allive-dead;
    
    countCHs=0;
    cluster=1;
%     CH selection
    for i=1:1:n
        if(S(i).E>=0)  %Checking threshold and nmbr of sleep nodes
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    S(i).type='C';
                    S(i).G=round(1/p)-1;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*(distance * distance ));
                    end
                end 
            end
        end
    end
    STATISTICS.COUNTCHS(r+1)=countCHs;
    %   Association of nodes
    for i=1:1:n
        if ( S(i).type=='N' && S(i).E>0)
            if(cluster-1>=1)
                min_dis=Inf;
                min_dis_cluster=0;
                for c=1:1:cluster-1
                    temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end
                min_dis;
                if (min_dis>do)
                    S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                    S(C(min_dis_cluster).id).E =S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                    packets_TO_CH=packets_TO_CH+1;
                
                    S(i).min_dis=min_dis;
                    S(i).min_dis_cluster=min_dis_cluster; 
     else
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
           if (min_dis>do)
               S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
           end
           if (min_dis<=do)
               S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
           end
           packets_TO_BS=packets_TO_BS+1;
            end
        end
    end    
                
    STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
    STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:n
    S1(i).xd=rand(1,1)*xm;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S1(i).yd=rand(1,1)*ym;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S1(i).G=0;                        % as the no of node that have been cluster head is zero 0
    S1(i).E=Eo;
    %initially there are no cluster heads only nodes
    S1(i).type='N';
     S(i).role=0;
end
S1(n+1).xd=sink.x;   %assume that base station is also a node sp total no of nodes is n and with base station  it is n+1
S1(n+1).yd=sink.y;
% % ETX=50*0.000000001;
% % ERX=50*0.0000000001;
ERX=50*0.0000000005;
 Efs=10*0.0000000000001;     %amplification energy when d is less than d0
Emp=0.0013*0.000000000000001;
countCHs1=0;         %the number of Stateflow objects in the current context.
cluster1=1;              %first cluster is selected
flag_first_dead1=0;
flag_teenth_dead1=0;
flag_all_dead1=0;
dead1=0;
first_dead1=0;
teenth_dead1=0;
all_dead1=0;
allive1=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS1=0;
packets_TO_CH1=0;
% counter for sleep nodes
s1=0;
th=0.00000000000001;
for r=0:1:rmax
    r
    if(mod(r, round(1/p) )==0) %remainder
        for i=1:1:n
            S1(i).G=0;            % it will assign to the nodes that have not been cluster head .
        end
    end
    
    dead1=0;
    for i=1:1:n
        
        if (S1(i).E<=0)
            dead1=dead1+1;
            
            if (dead1==1)
                if(flag_first_dead1==0)
                    first_dead1=r;
                    flag_first_dead1=1;
                end
            end
            
            if(dead1==0.1*n)
                if(flag_teenth_dead1==0)
                    teenth_dead1=r;
                    flag_teenth_dead1=1;
                end
            end
            if(dead1==n)
                if(flag_all_dead1==0)
                    all_dead1=r;
                    flag_all_dead1=1;
                end
            end
        end
%         check for sleep nodes
        if (S1(i).E<=th && S1(i).E>0)
            s1=s1+1;
        end
        
        if (S1(i).E>0)
            S1(i).type='N';
        end
    end
    STATISTICS.DEAD1(r+1)=dead1;
    STATISTICS.ALLIVE1(r+1)=allive1-dead1;
    
    countCHs1=0;
    cluster1=1;
%     CH selection
    for i=1:1:n
        if(S1(i).E>=th && s1<10)  %Checking threshold and nmbr of sleep nodes
            temp_rand=rand;
            if ( (S1(i).G)<=0)
                
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    countCHs1=countCHs1+1;
                    packets_TO_BS1=packets_TO_BS1+1;
                    PACKETS_TO_BS1(r+1)=packets_TO_BS1;
                    S1(i).type='C';
                    S1(i).G=round(1/p)-1;
                    C(cluster1).xd=S1(i).xd;
                    C(cluster1).yd=S1(i).yd;
                    distance=sqrt( (S1(i).xd-(S1(n+1).xd) )^2 + (S1(i).yd-(S1(n+1).yd) )^2 );
                     S(i).dts=distance;
                    C(cluster1).distance=distance;
                    C(cluster1).id=i;
                    X(cluster1)=S1(i).xd;
                    Y(cluster1)=S1(i).yd;
                    cluster1=cluster1+1;
                    distance;
                    if (distance>do)
                        S1(i).E=S1(i).E- ( (ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S1(i).E=S1(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*(distance * distance ));
                    end
                end 
            end
        end
    end
    STATISTICS.COUNTCHS1(r+1)=countCHs1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    X=3;
cluster2=2;
Th_cluster=0.0001;
B = sort(cluster1,'descend');
for i=1:X
    i=i+1;
    Larg=i;
   S1(i).role=max(Larg);
   for j=1:S1(i).role
       j=j+1;
      S1(i).dts=countCHs1;
       SCH=max(countCHs1);
   end
   B1 = sort(countCHs1,'descend');
   k=S1(i).role-Th_cluster;
   if k>0 
       for k=1:n
           B1;
           cluster2=cluster1+1;
       end
   end
   D=STATISTICS.COUNTCHS(r+1)-SCH;
    B1 = sort(D,'descend');
end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %   Association of nodes
    for i=1:1:n
        if ( S1(i).type=='N' && S1(i).E>=th && s1<10)
            if(cluster1-1>=1)
                min_dis=Inf;
                min_dis_cluster1=0;
                for c=1:1:cluster1-1
                    temp=min(min_dis,sqrt( (S1(i).xd-C(c).xd)^2 + (S1(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster1=c;
                    end
                end
                
                
                min_dis;
                if (min_dis>do)
                    S1(i).E=S1(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S1(i).E=S1(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                S1(C(min_dis_cluster1).id).E =S1(C(min_dis_cluster1).id).E- ( (ERX + EDA)*4000 );
                packets_TO_CH1=packets_TO_CH1+1;
                
                S1(i).min_dis=min_dis;
                S1(i).min_dis_cluster1=min_dis_cluster1;  
            end
        end
        %   Activation of sleep nodes when their nmbr exceed 9
        if(s1>=10 && S1(i).E>=0)
            if(cluster1-1>=1)
                min_dis=Inf;
                min_dis_cluster1=0;
                for c=1:1:cluster1-1
                    temp=min(min_dis,sqrt( (S1(i).xd-C(c).xd)^2 + (S1(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster1=c;
                    end
                end
                
                
                min_dis;
                if (min_dis>do)
                    S1(i).E=S1(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S1(i).E=S1(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                S1(C(min_dis_cluster1).id).E =S1(C(min_dis_cluster1).id).E- ( (ERX + EDA)*4000 );
                packets_TO_CH1=packets_TO_CH1+1;
                S1(i).min_dis=min_dis;
                S1(i).min_dis_cluster1=min_dis_cluster1;
            end
        end
        
    end
    STATISTICS.PACKETS_TO_CH1(r+1)=packets_TO_CH1;
    STATISTICS.PACKETS_TO_BS1(r+1)=packets_TO_BS1;
end
figure(1);
r=0:rmax;
plot(r,STATISTICS.ALLIVE(r+1),'-b', r,STATISTICS.ALLIVE1(r+1),'-r')
legend('LEACH   ','LEACH-USC');
xlabel('No. of Rounds (r)');
ylabel('No. of  Nodes Allive');

