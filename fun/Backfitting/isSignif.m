function [isSignif_NCs,timeSignif]=isSignif(results_ICA,perms,NCs,ind_0s)

for i=1:NCs
    sig=squeeze(results_ICA.signals(i,:,:)); 
    meansig(i,:)=mean(sig,2);
    up(i,:)=perms.thresholds.upper(i,:); 
    low(i,:)=perms.thresholds.lower(i,:);
    if((any(meansig(i,ind_0s:end)>up(i,ind_0s:end)))||(any(meansig(i,ind_0s:end)<low(i,ind_0s:end))))
        isSignif_NCs{i}=true;
        upSignif=find(meansig(i,ind_0s:end)>up(i,ind_0s:end));
        lowSignif=find(meansig(i,ind_0s:end)<low(i,ind_0s:end));
        indSignif{i}=[upSignif+ind_0s-1 lowSignif+ind_0s-1];
        timeSignif{i}=results_ICA.time(indSignif{i});
        timeSignif{i}=sort(timeSignif{i},'ascend');
    else
        isSignif_NCs{i}=false;
        timeSignif{i}=[];
    end
    
end
