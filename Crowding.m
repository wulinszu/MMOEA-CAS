function [crowd rank]=Crowding(Pop)
    [N D]=size(Pop);
    K=N-1;
    Z = min(Pop,[],1);
    Zmax = max(Pop,[],1);
    pop=(Pop-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance=pdist2(pop,pop);
    [value,index]=sort(distance,2);
    crowd=K./sum(1./value(:,2:K),2);
    [~,rank]=sort(crowd,'ascend');
end