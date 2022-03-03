function results = go_decomposeConnectome_SS_ICA_noload(opt,varagin)

% This code takes as input opt structure that should contains the basic info for source separation:

% opt.data: cell of length nb of subject, containing cmat structure output from connectivity code (go_dynamicFC)
% opt.n_parcels: number of parcels or ROIs used
% opt.NCs: number of Results Components
% opt.val: int to index the method of ICA to be used (1:ICA-JADE, 2:ICA-InfoMax, 3:ICA-SOBI, 4:ICA-FastICA, 5:ICA-CoM2, 6:ICA-PSAUD)
% opt.cmat_all: concatenated cmat data for all subj and trials;
% opt.cmat_time: time vector;
% opt.index: index of ROIs reduced;
% opt.sub_trials: nb of trials for each subj;

% results: output structure containing results of Source Separation ICA Decomposition (connectivities maps and temporal signals)

% This code was originally developped by Judie Tabbal.
% contact: judytabal95@gmail.com

opt.NCs            = ft_getopt(opt, 'NCs', 10);
opt.approach        = ft_getopt(opt, 'approach', 'defl');

if ~isfield(opt,'data');
    error('List of connectime .mat files not present!');
end

% add the fastica toolbox to path if its not already present
tmp = which('fastica');
if isempty(tmp)
   [status] = ft_hastoolbox('fastica', 1, 0);
   if ~status
       error('fastica toolbox failed to load correctly')
   end
end

% load and concatenate all the data for ICA. 
n_subs      = length(opt.data);
cmat_all    = [];

disp('Loading and concatenating connectomes :')
ft_progress('init', 'text', 'Please wait...')
pause(0.01)

if(opt.val==1) %ICA-JADE
    [B] = jader(opt.cmat_all, opt.NCs);
    A=B'; 
    sig=B*opt.cmat_all; 
elseif(opt.val==2) %ICA-InfoMax
    [weights,sphere] = runica(opt.cmat_all,'pca',opt.NCs,'ncomps',opt.NCs);
    B=weights*sphere;
    A=B'; 
    sig=B*opt.cmat_all; 
elseif(opt.val==3) %ICA-SOBI
    [A, sig] = sobi_edit(opt.cmat_all, opt.NCs);
elseif(opt.val==4) %ICA-FastICA
%     [sig, A, W] = fastica(opt.cmat_all, 'g', 'tanh','lastEig', opt.NPCs, 'numOfIC', opt.NCs,'approach','defl');
    [iq, A, W, sig, sR]=icasso(opt.cmat_all,100,'g','tanh','lastEig', opt.NCs,'numOfIC', opt.NCs,'approach','defl','vis','off');
elseif(opt.val==5) %ICA-CoM2
    [F, delta] = aci(opt.cmat_all,opt.NCs);
    A=F; 
    sig=F'*opt.cmat_all; 
elseif(opt.val==6) %ICA-PSAUD
    tau=1;
    P=opt.NCs;
    Nit=20;
    alpha_min=0;
    alpha_max=50;
    [F,S]=P_SAUD(opt.cmat_all,tau,P,Nit,alpha_min,alpha_max);
    A=F; %estimated mixing matrix
    sig=S; %estimated signal matrix
end

% post processing of the ICA results 
tmp = reshape(sig,opt.NCs,size(sig,2)/sum(opt.sub_trials),sum(opt.sub_trials));
results.signals = tmp;

for ii = 1:opt.NCs
    tmp = zeros(opt.n_parcels,opt.n_parcels);
    tmp(opt.index) = A(:,ii);
    tmp = tmp + transpose(tmp);
    results.maps(:,:,ii) = tmp;
end

results.NCs        = opt.NCs;
results.n_trials    = sum(opt.sub_trials);
results.sub_trials  = opt.sub_trials;
results.time        = round(opt.cmat_time*100)/100;
results.sig=sig;
results.mixing=A;
results.cmat_all=opt.cmat_all;