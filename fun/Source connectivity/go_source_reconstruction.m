function filters=go_source_reconstruction(opt,data,subvol,subgrid,sources_orients)

%Calculate Lead fields (Forward Model)
cfg             = [];
cfg.headmodel   = subvol;
cfg.elec        = data.elec;
cfg.grid.pos    = subgrid.pos;
cfg.normalise   = 'yes';
cfg.rankreduce  = 3; 
lf              = ft_prepare_leadfield(cfg); 

   
%Timelock for noise covariance estimation
cfg                     = [];
cfg.covariance          = 'yes';
cfg.window              = [-opt.prestim 0]; %baseline for noise cov (computed from prestim in seconds to 0)
tlk_noise                = ft_timelockanalysis(cfg,data);
noise_cov=tlk_noise.cov;

filters = ComputeWMNE(noise_cov,cell2mat(lf.leadfield),lf.pos,sources_orients,opt.weightExp,opt.weightLimit,opt.SNR);

% % uncomment below only for source plot visualization (indicate interval of time to average in cfg.latency)
% cfg                 = [];
% cfg.method          = 'mne';
% cfg.grid            = lf;
% cfg.elec            = data.elec;
% cfg.headmodel       = subvol;
% cfg.mne.snr         = 3; %regularisation param for noise cov matrix 
% src                 = ft_sourceanalysis(cfg,data);
%
% nROIs=68;
% nSamples=1168;
% nTrials=200;
% scout_source=zeros(nROIs,nSamples);
% for tr=1:nTrials
%     scout_source=scout_source+abs(filters*data.trial{tr}); 
% end
% scout_source=scout_source./nTrials;
% src.aa=scout_source;
% 
% cfg = [];
% cfg.funparameter = 'aa';
% cfg.maskparameter = 'aa';
% cfg.method = 'surface';
% cfg.surffile= 'D:\Brainstorm\brainstorm_db\Protocol05\anat\Subject01\tess_cortex_pial_low.mat';
% cfg.latency = [0.190 0.340];
% cfg.avgovertime = 'yes';
% max_avg=max(mean(src.aa(:,floor((cfg.latency(1)+0.3)*data.fsample):floor((cfg.latency(2)+0.3)*data.fsample)),2));
% cfg.opacitylim = [0.2*max_avg max_avg];
% ft_sourceplot(cfg, src);