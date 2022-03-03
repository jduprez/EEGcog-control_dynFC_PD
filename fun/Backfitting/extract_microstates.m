function [microparams_grp]=extract_microstates(opt,corr_tw,ind_tw,cmat_allstr,n_subs)

% Inputs:
%   corr_tw     : correlation values between spatial maps of each subject and at each temporal window with each of the significant maps states of both groups (output of do_backfitting.m)
%   ind_tw      : the index of significant state with the maximum correlation value (for each subject and at each temporal window) (output of do_backfitting.m)
%   cmat_allstr : output cell of length (n_subs) containing loaded cmat results for each subject

%   n_subs      : number of subjects 
%   opt         : structure containing the following field:
%                 opt.cCAT      : number of significant kept states 
%                 opt.ind_0s    : onset time (in sec)
%                 opt.totaltime : total time duration after onset (in sec)
%                 opt.deltat    : difference time between two time windows (in sec)
%
% Output:
%   microparams_grp : output structure of microstats parameters relative to the group, with the following fields:
%                     microparams_grp.fraction_covtime : fraction of time a given state is active
%                     microparams_grp.freq_occurence   : average number of times per second a state is dominant
%                     microparams_grp.avg_lifespan     : average duration of a given state 
%                     microparams_grp.GEV              : global explained variance of a given state
%                     microparams_grp.TR               : transition probabilities between states
%                     microparams_grp.TRsym0           : transition probabilities between states (symmetric version with diagonal = 0)



% loop over all significant kept states
for istate=1:opt.cCAT
    for isub=1:n_subs
        vector_corr=corr_tw{isub}(istate,opt.ind_0s:end);
        vector=ind_tw{isub}(opt.ind_0s:end);
        ind_istate=find(vector==istate);
        if(isempty(ind_istate))
            fraction_covtime{istate}(isub,1)=0;
            freq_occurence{istate}(isub,1)=0;
            avg_lifespan{istate}(isub,1)=0;
            GEV{istate}(isub,1)=0;
        else
            %calcul fraction coverage time
            fraction_covtime{istate}(isub,1)=length(ind_istate)/length(vector);
            
            %calcul freq occurence
            diff_ind_istate=diff(ind_istate)==1;
            ind0_diff_ind_istate=find(diff_ind_istate==0);
            freq_count=length(ind0_diff_ind_istate)+1;
            freq_occurence{istate}(isub,1)=freq_count/opt.totaltime;
            
            %calcul avg lifespan
            new_vec_ind=[0 diff_ind_istate];
            fc=1;
            dur=ones(1,freq_count);
            for ii=2:length(new_vec_ind)
                if(new_vec_ind(ii)==0)
                    fc=fc+1;        
                elseif(new_vec_ind(ii)==1)
                    dur(fc)=dur(fc)+1;
                end
            end
            avg_lifespan{istate}(isub,1)=mean(dur)*opt.deltat;
            
            %calcul GEV (Global Explained Variance)
            cmat=cmat_allstr{isub};
            nROI=size(cmat.connectivity{1},1);
            for itr=1:length(cmat.connectivity)
                conn_alltr(:,:,:,itr)=cmat.connectivity{itr};
            end
            conn_avgtr=mean(conn_alltr,4);
            net=conn_avgtr(:,:,opt.ind_0s:end);
            GEV_num=var(var(reshape(net(:,:,vector==istate),nROI*nROI,sum(vector==istate)).*vector_corr(vector==istate)));
            GEV_denum=var(var(reshape(net(:,:,vector==istate),nROI*nROI,sum(vector==istate))));
            GEV{istate}(isub,1)=GEV_num./GEV_denum;
        end
    end
    
%calcul TransitionProba between states
for isub=1:n_subs %n_controls
    TR_i=zeros(opt.cCAT);
    vector=ind_tw{isub}(opt.ind_0s:end);
    for ii=1:length(vector)-1
        st1=vector(ii);st2=vector(ii+1);
        TR_i(st1,st2)=TR_i(st1,st2)+1;
    end
    TR_i_rowsum=sum(TR_i,2);
    TR_i_rowsum(TR_i_rowsum==0)=-inf;
    TR_i_freq=TR_i./repmat(TR_i_rowsum,1,opt.cCAT);
    TR{isub}=TR_i_freq;
    %TR symmetric with diag=0: TRsym0
    TRsym0_i=(TR_i+TR_i')-diag(diag(TR_i+TR_i'));
    TRsym0_i_rowsum=sum(TRsym0_i,2);
    TRsym0_i_rowsum(TRsym0_i_rowsum==0)=-inf;
    TRsym0_i_freq=TRsym0_i./repmat(TRsym0_i_rowsum,1,opt.cCAT);
    TRsym0{isub}=TRsym0_i_freq;
end

microparams_grp.fraction_covtime=fraction_covtime;
microparams_grp.freq_occurence=freq_occurence;
microparams_grp.avg_lifespan=avg_lifespan;
microparams_grp.GEV=GEV;
microparams_grp.TR=TR;
microparams_grp.TRsym0=TRsym0;

end

