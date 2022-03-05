function [isSignif_NCs,timeSignif]=isSignif_cycle(results_ICA,perms,NCs,ind_0s,d_cy)

% This function extract significant states with corresponding significance time 
% (based on null distribution and 3+ significant cycles)
% 
% Inputs: 
%   results_ICA : ICA results structure (output of go_decomposeConnectome_SS_ICA.m)
%   perms       : permutations results (output of go_testNetworks_general.m)
%   NCs         : number of all extracted independent components
%   ind_0s      : index sample of onset time
%   d_cy        : duration of 3 cycles 
%
% Outputs:
%   isSignif_NCs : boolean varibale (true/false) that indicate for each state if it is significant to keep (true) or not (false)
%   timeSignif   : array containing the times at which the state is significant



for i=1:NCs
    
    % 1. Significance based on null distribution 
    sig = squeeze(results_ICA.signals(i,:,:));
    meansig(i,:) = mean(sig,2);
    up(i,:) = perms.thresholds.upper(i,:);
    low(i,:) = perms.thresholds.lower(i,:);
    
    if((any(meansig(i,ind_0s:end)>up(i,ind_0s:end)))||(any(meansig(i,ind_0s:end)<low(i,ind_0s:end))))
        isSignif_NCs{i} = true;
        upSignif = find(meansig(i,ind_0s:end)>up(i,ind_0s:end));
        lowSignif = find(meansig(i,ind_0s:end)<low(i,ind_0s:end));
        indSignif{i} = [upSignif+ind_0s-1 lowSignif+ind_0s-1];
        timeSignif{i} = results_ICA.time(indSignif{i});
        timeSignif{i} = sort(timeSignif{i},'ascend');
    else
        isSignif_NCs{i} = false;
        timeSignif{i} = [];
    end
    
    % 2. Significance based on number of significant cycles 
    
    % indexes of the  sig values
    [r, idx] = find(meansig(i,ind_0s:end)>up(i, ind_0s:end) | meansig(i,ind_0s:end)<low(i, ind_0s:end));
    
    time = results_ICA.time(ind_0s:end);
    
    D = diff([0, diff(idx)==1, 0]);
    first = idx(D>0);
    last = idx(D<0);
    
    % there is 109 samples from -0.62 s to 1.1 s
    
    decision = zeros(1, length(first));
    
    for wdwi = 1:length(first)
        if ((time(last(wdwi)) - time(first(wdwi)))*1000)>= d_cy
            decision(wdwi) = 1;
        else
            decision(wdwi) = 0;
        end
    end
    
    if sum(decision)>=1
        isSignif_NCs{i} = true;
    else
        isSignif_NCs{i}=false;
    end
    
    
end
