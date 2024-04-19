close all;
clear all;
%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%
%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
% x and y coordinates of the sink
sink.x=0.5*xm;
sink.y=0.5*ym;
%Number of nodes in the field
n=100;
%Optimal election probability of a node
%to be come cluster head
p=0.1;
%Energy model (ALL VALUES IN JOULES)
%Intial Energy
Eo=0.5;
%Eelec=Etx=Erc
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit amplifier types
Efs=10*0.00000000001;
Emp=0.0013*0.00000000001;
%Data aggreagation energy
EDA=5*0.000000001;
%%Values for Hetereogeneity
%Percentage of nodes than are advanced 
m=1.5;
%\alpha
a=0;
Ave_CH2=0;
sum2=0;
count_ch2=0;
Throughput2=0;
x=0.3;
b=0;
rmax=5000;
Ave_CHs=0;
sums=0;
count_chs=0;
Throughputs=0;
Packet=4000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PARAMETERS%%%%%%%%%%%%%%%%%

%Computaion of do
do=sqrt(Efs/Emp);

%Creation of the random sensor network
% figure(1);


rand('seed',13)
rng;

 for i2=1:1:n
    S2(i2).xd=rand(1,1)*xm;
    XR(i2)=S2(i2).xd;
    S2(i2).yd=rand(1,1)*ym;
    YR(i2)=S2(i2).yd;
    S2(i2).G=0;
    S2(i2).E=0;
    %intially there are no cluster heads only nodes
    S2(i2).type='N';
    S2(i2).id=i2;
    keep(i2)=i2;
    temp_rnd0=i2; 
     %Random election of normal nodes
    if(temp_rnd0>=(x+m)*n+1)
        S2(i2).E=Eo;
        S2(i2).ENERGY=0;
%         plot(SS(I).xd,S(i).yd,'0');
%         hold on;
    end
    
    %Random election of intermediate  nodes
    if(temp_rnd0<(m+x)*n+1)&&(temp_rnd0>m*n)
        S2(i2).E=Eo*(1+b);
        S2(i2).ENERGY=0.5;
%         plot(S(i).xd,S(i).yd,'+');
%         hold on;
    end
      %Random election of advanced nodes
       if(temp_rnd0<m*n+1)
        S2(i2).E=Eo*(1+a);
        S2(i2).ENERGY=1;
       end
    
end

S2(n+1).xd=sink.x;
S2(n+1).yd=sink.y;

S2(n+1).xd=sink.x;
S2(n+1).yd=sink.y;
% plot (S(n+1).xd,S(n+1).yd,'x');

%First Iteration


%counter for CHs
countCHs2=0;
%counter for CHs per round 
rcountCHs2=0;
cluster2=1;

countCHs2;
rcountCHs2=rcountCHs2+countCHs2;
flag_first_dead2=0;


for r2=0:1:rmax
    r2

  %Operation for epoch
    if(mod(r2, round(1/p) )==0)
        for i2=1:1:n
            S2(i2).G=0;
            S2(i2).c1=0;
        end
    end
    
%  hold off;

 %Number of dead nodes
 dead2=0;
 %Number of dead Advanced Nodes
 dead_a2=0;
 %Number of dead Normal Nodes
 dead_n2=0;
 dead_i2=0;
 %counter for bit transmitted to base station and to cluster head
 packets_TO_BS2=0;
 packets_TO_CH2=0;
 
 %counter for bit transmitted to base station and to cluster head
 %per round
 PACKETS_TO_CH2(r2+1)=0;
 PACKETS_TO_BS2(r2+1)=0;
 
 %figure(1);     
    
    
    
    for i2=1:1:n
     total_energy(i2)=S2(i2).E;
     %checking if there is a dead node
     if(S2(i2).E<=0)
         %plot(S(i).xd,S(i).yd,'red..');
         dead2=dead2+1;
         if(S2(i2).ENERGY==1)
               dead_a2=dead_a2+1;
         end
         if(S2(i2).ENERGY==0)
               dead_n2=dead_n2 +1;
         end
%          hold on;
     end
if S2(i2).E>0
         S2(i2).type='N';
         if (S2(i2).ENERGY==0)
         %plot(S(i).xd,S(i).yd,'0');
         end
          if (S2(i2).ENERGY==0.5)
               %plot(S(i).xd,S(i).yd,'p');
          end
         if (S2(i2).ENERGY==1)
         %plot(S(i).xd,S(i).yd,'+');
         end
end
     end

Total_energy=sum(total_energy);
Totalnetwor_energy(r2+1)=Total_energy;
 STATISTICS(r2+1).DEAD = dead2;
 DEAD2(r2+1) = dead2;
 DEAD_N2(r2+1) = dead_n2;
 DEAD_i2(r2+1) = dead_i2;
 DEAD_A2(r2+1) = dead_a2;
%When the first node dies
 if(dead2==1)
     if(flag_first_dead2==0)
         first_dead2=r2;
         flag_first_dead2=1;
     end
 end

%no of alive nodes
alive2=0;
%no of alive  normal nodes
alive_n2=0;
%no of alive advance nodes
alive_a2=0;
alive_i2=0;

for i2=1:1:n
%checking number of alive node per round
   if(S2(i2).E>0)
         alive2=alive2+1;
     if (S2(i2).ENERGY==1)
              alive_a2=alive_a2+1;
     end
        if (S2(i2).ENERGY==0)
              alive_n2=alive_n2+1;
        end
           if (S2(i2).ENERGY==0.5)
                alive_i2=alive_i2 +1;
           end
   end
%checking nodes status
      if(S2(i2).E>0)
           nodes_status2=1;
      end
        if(S2(i2).E<0)
            nodes_status2=0;
        end
%    end
STATISTICS(i2).Status=nodes_status2;
Status2(i2)=nodes_status2;
ASTATISTICTS2(r2+1).Live=alive2;
Live2(r2+1)=alive2;
Live_n2(r2+1)=alive_n2;
Live_i2(r2+1)=alive_i2;
Live_a2(r2+1)=alive_a2;
end
%checkingfor last dead or alive node
for i2=1:1:n
    if(alive2==1&&S2(i2).E>0)
        if(S2(i2).ENERGY==1||S2(i2).ENERGY==0||S2(i2).ENERGY==0.5)
  last_dead2=r2;
instability2=last_dead2-first_dead2;
flag_last_dead2=1;
        end
    end
end
 countCHs2=0;
 cluster2=1;
 for ii2=1:1:n
     if(S2(ii2).E>0)
         temp_rand2=rand;
         if ( (S2(ii2).G)<=0)
             
             %Election of cluster heads
                distance2=sqrt( (S2(ii2).xd-(S2(n+1).xd) )^2+ (S2(ii2).yd-(S2(n+1).yd) )^2 );
             if(temp_rand2<=(p/(1-p*mod(r2,round(1/p)))+(sqrt(2)- distance2/(S2(ii2).xd+S2(ii2).yd)/4)))
                  packets_TO_BS2=packets_TO_BS2+1;  
                 PACKETS_TO_BS2(r2+1)=packets_TO_BS2;
                 
                 S2(ii2).type='C';
                 S2(ii2).G=round(1/p)-1;
                 C2(cluster2).xd=S2(ii2).xd;
                 C2(cluster2).yd=S2(ii2).yd;
           
                  
                  distance2=sqrt( (S2(ii2).xd-(S2(n+1).xd) )^2+ (S2(ii2).yd-(S2(n+1).yd) )^2 );
                  C2(cluster2).distance=distance2;
                  C2(cluster2).id=S2(ii2).id;
                  X2(cluster2)=S2(ii2).xd;
                  Y2(cluster2)=S2(ii2).yd;
                  cluster2=cluster2+1;
                  %calculation of Energy dissipated
                  distance2;
             end
         end
     end
 end
 %checking average number of Cluster HEads per epoch
                  sum2=sum2+(cluster2-1);
                  count_ch2=count_ch2+1;
                  I2=100;
                  if count_ch2==10
                      Ave_CH2=(sum2*0.1)/(1+(m*a));
                      Throughput2= Ave_CH2*4;
                      STATSTISTICS(r2+1).ave_clustHd=Ave_CH2;
                      ave_ch2(r2+1)=Ave_CH2;
                      STATSTISTICS(r2+1).throughput=Throughput2;
                      Clust_throughput(r2+1)=Throughput2;
                      if Ave_CH2~=0
                          count_object(I2)=Ave_CH2;
                      end
                      Ave_CH2=0;
                      sum2=0;
                      count_ch2=0;
                  end
               STATISTICS(r2+1).CLUSTERHEADS=cluster2-1;
               CLUSTERHS2(r2+1)=cluster2-1;
               %election of cluster head for normal nodes
               for i2=1:1:n
                    if (S2(i2).type == 'N' && S2(i2).E>0)
                       if(cluster2-1>=1)
                           min_dis2=9999;
                           min_dis_cluster2=1;
                           for c=1:1:cluster2-1
                               temp2=min(min_dis2,sqrt((S2(i2).xd-C2(c).xd)^2+(S2(i2).yd-C2(c).yd)^2));
                               if(temp2<min_dis2)
                                   min_dis2=temp2;
                                   min_dis_cluster2=c;
                               end
                           end
           %energy dissipated by associated Vice Cluster Head(subcluster finding)
           min_dis2;
            if (min_dis2<=do)
               S2(i2).E=S2(i2).E- (ETX*(4000) + Efs*4000*( min_dis2 * min_dis2));
           end
         if (min_dis2>do) 
             for i=1:1:c-1
                 %-------search crieteria:
                 temps2=sqrt((S2(i).xd-C2(c).xd)^2+(S2(i).yd-C2(c).yd)^2);
               min_dis2=9999999;
               if temps2>do
                   for j=1:1:c-1
                       find_near_subcluster=min(min_dis2,sqrt((S2(i).xd-S2(j).xd)^2+(S2(i).yd-S2(j).yd)^2));
                       if(find_near_subcluster< min_dis2)
                           min_dis2=find_near_subcluster;
                       end
                   end
                 S2(i).E=S2(i).E- (ETX*(4000) + Efs*4000*( min_dis2 * min_dis2));  
                      
               end
               
             end
         end
         %ENERGY Dissipated
         if(min_dis2>0)
             S2(C2(min_dis_cluster2).id).E=S2(C2(min_dis_cluster2).id).E-(ERX+EDA)*4000;
             PACKETS_TO_CH2(r2+1)=n-dead2-cluster2+1;
         end
         S2(i2).min_dis=min_dis2;
         S2(i2).min_dis_cluster=min_dis_cluster2 ;
                       end
                    end
               end
    %election of next hop Dynamic election of Vice Cluster head
    countCHs2;
    rcountCHs2=rcountCHs2+countCHs2;
    for c=1:1:cluster2-1        
      tempp2=sqrt((S2(n+1).xd-C2(c).xd)^2+(S2(n+1).yd-C2(c).yd)^2);              
        diss(c)=tempp2;
        C2(c).id;
        if tempp2<do
            S2(C2(c).id).E=S2(C2(c).id).E-((ETX+EDA)*(4000)+Efs*4000*(tempp2*tempp2));
        end
        All_dis(c)=tempp2;
    end
    [diss,idd]=sort(diss,'ascend');
    for i=1:1:cluster2-1
        temps2=(sqrt((S2(n+1).xd-C2(i).xd)^2+(S2(n+1).yd-C2(i).yd)^2));
        min_dis2=9999999;
        if temps2>do
            for j=1:1:cluster2-1
                       find_near_cluster=min(min_dis2,sqrt((C2(i).xd-C2(j).xd)^2+(C2(i).yd-C2(j).yd)^2));
                       if(find_near_subcluster< min_dis2)
                           min_dis2=find_near_cluster;
                       end
            end 
            S2(C2(i).id).E=S2(C2(i).id).E-(ETX*(4000)+Efs*4000*(min_dis2*min_dis2));
        end
    end
     
 
end     
figure(1)
plot(1:5001,Live2)
xlabel('No. of Rounds');
ylabel('No. of Nodes Alive')
title('Alive Nodes');
figure(2) 
plot(sort(Totalnetwor_energy/100,'ascend'),'rd-');
xlabel('rounds');
ylabel('consumed energy');
title('Energy consumption for nodes');
    
    
            
            
            
            