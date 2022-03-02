function [cmat_all,index,sub_trials,cmat_time] = load_and_concatenate(opt,varagin)

% This code takes as input opt structure that should contains the basic info for source separation:

% opt.data: cell of length nb of subject, containing cmat structure output from connectivity code (go_dynamicFC)
% opt.n_parcels: number of parcels or ROIs used
% opt.NCs: number of Results Components

% load and concatenate all the data for ICA. 
n_subs      = length(opt.data);
cmat_all    = [];

disp('Loading and concatenating connectomes :')
ft_progress('init', 'text', 'Please wait...')
pause(0.01)

for ii = 1:n_subs
    
    ft_progress(ii/n_subs, 'Processing connectome %d of %d', ii, n_subs);   
    load(opt.data{ii});
    
    cmat_sub = [];
    cmat_2d  = [];
    cmat_red = [];
    
    % loop across trials
    for jj = 1:cmat.n_trials
        tmp = cmat.connectivity{jj};
        cmat_sub = cat(3,cmat_sub,tmp);
    end   
    
    % Do a DC correction on each connection
    cmat_sub = cmat_sub-repmat(mean(cmat_sub,3),[1 1 size(cmat_sub,3)]);   
    % flatten tensor into 2.5d (as the machine learning kids call it)
    cmat_2d = reshape(cmat_sub,opt.n_parcels*opt.n_parcels,length(cmat.time)*cmat.n_trials);
    % as each slice of the tensor is symmetric, remove redundant data to save some memory
    index = find(tril(squeeze(cmat_sub(:,:,1)))~=0);
    cmat_red = cmat_2d(index,:);     
    % concatenate onto the group level connectome data
    cmat_all = cat(2,cmat_all,cmat_red);
    
    sub_trials(ii) = cmat.n_trials; 
    
    cmat_time=cmat.time;
end

disp('DONE')