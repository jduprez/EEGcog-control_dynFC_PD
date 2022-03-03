function [corr_tw,max_tw,ind_tw,cmat_allstr]=do_backfitting(opt,cmat_list,n_subs)

% This function assign the most similar spatial map among significant
% states (extracted at the group-level) to each of subject-level spatial
% maps and at each temporal window
%
% Inputs:
%   cmat_list : cell of length (nb_subs) containing the path of saved connectivity PLV
%   n_subs    : number of subjects
%   opt       : input structure with the following fields:
%               opt.threshnet_meth  : apply threshold on maps before similarity calculation (only 'no' option is implemented in the code
%               opt.corr_meth       : apply correlation as spatial similarity measure (only 'corr2' option is implemented in the code)
%               opt.cCAT             : number of significant kept states 
%               opt.states_maps_all : combined HC followed by PD significant maps, 3D matrix of dim [nROIs,nROIs,cHC+cPD];
%
% Outputs: 
%   corr_tw     : correlation values between spatial maps of each subject and at each temporal window with each of the significant maps states of both groups
%   max_tw      : the maximum correlation value among significant states (for each subject and at each temporal window)
%   ind_tw      : the index of significant state with the maximum correlation value (for each subject and at each temporal window)
%   cmat_allstr : output cell of length (n_subs) containing loaded cmat results for each subject



% loop over subjects
for isub=1:n_subs
    load(cmat_list{isub});
    cmat.connectivity = cellfun(@single, cmat.connectivity, 'un', 0); % put to single to save space
    cmat_allstr{isub}=cmat;
    nROI=size(cmat.connectivity{1},1);
    for itr=1:length(cmat.connectivity)
        conn_alltr(:,:,:,itr)=cmat.connectivity{itr};
    end
    clearvars cmat
    conn_avgtr=mean(double(conn_alltr),4); % put to double again just during computation of the mean !
    n_windows=size(conn_avgtr,3);
    for itw=1:n_windows
        for istate=1:opt.cCAT
            net1_comp=conn_avgtr(:,:,itw);
            net2_comp=opt.states_maps_all(:,:,istate);
            corr_tw{isub}(istate,itw)=corr2(net1_comp,net2_comp);
        end
    end
    for itw=1:n_windows
        [max_tw{isub}(itw),ind_tw{isub}(itw)]=max(corr_tw{isub}(:,itw));
    end
    clear cmat
end


end

