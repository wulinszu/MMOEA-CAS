function NSGAII7_1(Global)
% <algorithm> <H-N>
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool=randperm(Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        %联合种群
        union=[Population Offspring];
        %global_archive保存全局最优解
        [FrontNo,~]=NDSort(union.objs,length(union));
        global_A = find(FrontNo==1);
        %基于决策空间邻域关系的聚类
        range=max(union.decs,[],1)-min(union.decs,[],1);
        epsilon=mean(range)*0.1;
        K=3;
        IDX=DBSCAN(union.decs,epsilon, K);
        local_A=[];
        for i=1:max(IDX)
            cluster=find(IDX==i);
            [FrontNo1,~]=NDSort(union(cluster).objs,length(cluster));
            %local_archive存储各个簇中的非支配解
            local_A=[local_A; cluster(FrontNo1==1)];
        end
        temp1=setdiff(global_A,local_A)';
        Q=[temp1; local_A];
        S=randperm(2*Global.N);
        U=setdiff(S,Q);
        temp=setdiff(local_A,global_A)';
%         rate=abs(length(global_A)-length(local_A))/length(local_A);
        C=[];
        i=1;
        [W]=UniformPoint(Global.N,Global.M);
        while  length(Q)<Global.N
            i=i+1;
            s=U((FrontNo(U)==i))';
            Q=[Q; s]; 
        end
        newpop=union(Q);
        z=min(newpop.objs,[],1);
        Z=max(newpop.objs,[],1);
        [~,Region] = max(1-pdist2((newpop.objs-z)./(Z-z),W,'cosine'),[],2);
        for i=1:length(newpop)-Global.N
            [crowd,rank]=Crowding(newpop.decs);
            c=hist(Region,1:size(W,1));
            %a中的索引表示最拥挤的子区域
            a=find(c==max(c));
            s=[];
            for j=1:length(a)
                s=[s;find(Region==a(j))];
            end
            s_old=s;
            ccc=0;
            while length(s)>0
                kk=IDX(s(1));
                sss=find(IDX(s)==kk);
                if length(sss)>ccc
                   ccc= length(sss);
                   index=kk;
                end
                s=setdiff(s,s(sss));
            end
            sss=find(IDX(s_old)==index);
            [~,s1]=sort(crowd(s_old(sss)),'ascend');
            x=s_old(sss(s1(1)));
            newpop(x)=[];
            Region(x)= [];
        end
        Population=newpop;
    end
end

