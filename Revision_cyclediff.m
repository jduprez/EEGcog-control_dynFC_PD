%% Info TBD

%% 1.1. Define basic param and add basic folders to matlab path


path_base = 'F:\WJD\Simon Dynamic FC';
mri_template = 'icbm'; % 'colin' for Colin27 or 'icbm' for ICBM152
atlas = 'destrieux'; % 'desikan' for Desikan68 atlas or 'destrieux' for Destrieux148 atlas
srate = 1000; % Define sampling rate
frequency = 'beta';

% Automatically add folders and subfolders (Code and Inputs) necessary for
% loading variables and executing codes
addpath(genpath(['F:\WJD\Simon Dynamic FC\input4code\' mri_template]));
addpath(genpath('C:\GitHub\EEGcog-control_dynFC_PD'));

% Add openmeeg path
setenv('PATH', [path_base '\Toolboxes\OpenMEEG\bin'])
addpath([path_base '\Toolboxes\fieldtrip-20190224\external\openmeeg'])

% Add default fieldtrip path
tmp = which('ft_defaults');
if isempty(tmp)
    addpath([path_base '\Toolboxes\fieldtrip-20190224']); %add defaults fieldtrip
    ft_defaults
end



%% 1.2. Define wMNE specification (prefer to keep these default params values)

source.weightExp = 0.5; % param for depth weighting
source.weightLimit = 10; % param for depth weighting limit
source.SNR = 3; % signal to Noise Ratio
% Load already computed and saved variables (as we are dealing with template
% MRI; these variables were computed for once as in code_for_inputs.m function
% and loaded here to be used for all subjects

load([path_base '\input4code\' mri_template '\elec_BS_' mri_template '_199ch.mat']); % fieldtrip format electrode already computed from modified EGI257 system: extracted from BrainStorm
load([path_base '\input4code\' mri_template '\' atlas '\' atlas '_scs_' mri_template '.mat']); % scout structure to extract sources centroids positions orientations (in SCS or CTF coord): calculated in Section 2 of 'code_for_inputs'
load([path_base '\input4code\' mri_template '\' atlas '\subgrid.mat']); % fieldtrip format of the grid specific to the mri and scout used (refer to code_for_inputs.m for details about computation): calculated in Section 3 of 'code_for_inputs'
load([path_base '\input4code\' atlas '_labels.mat']); % cell labels for the scout regions: extracted from BrainStorm


%% 1.3. Define dynamic Functional Connectivity (dFC) specification
for ncycle = [4 6 8]
    for perc = [70 80 90]
        dFC.conn_method = 'plv_dyn'; %'plv_dyn' for windowedPLV or 'wPLI' for windowedwPLI or 'plv_inst_pn' for instantaneous PLV
        if strcmp('gamma', frequency)
            dFC.bpfreq = [30 45]; %frequency band of interest
            central_freq = dFC.bpfreq(1)+(dFC.bpfreq(2)-dFC.bpfreq(1))/2;
            dFC.window.size = ncycle/central_freq; % sliding window length in seconds (for example calculated for 6cycles,CentralFreq=35 ==> 6/35=0.17s), in case of 'plv_inst_pn' this input is meaningless
            dFC.window.step = (100-perc)/100*dFC.window.size; % step between windows in seconds (for example 90% overlapping=10/100*window_size), in case of 'plv_inst_pn' this input is meaningless
            disp(['window size is ' num2str(dFC.window.size) ' and step is ' num2str(dFC.window.step)])
            
        elseif strcmp('beta', frequency)
            dFC.bpfreq = [12 25]; %frequency band of interest
            central_freq = dFC.bpfreq(1)+(dFC.bpfreq(2)-dFC.bpfreq(1))/2;
            dFC.window.size = ncycle/central_freq; % sliding window length in seconds (for example calculated for 6cycles,CentralFreq=35 ==> 6/35=0.17s), in case of 'plv_inst_pn' this input is meaningless
            dFC.window.step = (100-perc)/100*dFC.window.size; % step between windows in seconds (for example 90% overlapping=10/100*window_size), in case of 'plv_inst_pn' this input is meaningless
            disp(['window size is ' num2str(dFC.window.size) ' and step is ' num2str(dFC.window.step)])
            
        else
            error('Frequency is either absent or incompatible')
        end
        %% 3. Main code to compute sources+connectivities
        
        cd('F:\WJD\Simon Dynamic FC\Data')
        ftrialname = dir('*_incong.mat');
        
        %%%%%%%%%%%%%%% Manipulate these parameters appropriatly %%%%%%%%%%%%%%%
        
        nb_sub=length(ftrialname); %number of subjects
        fs=1000; %sampling frequency
        pre_samples=700; %nb samples before trial onset
        post_samples=1200; %nb samples after trial onset
        nb_samples=1901; %total number of samples in each trial
        nb_chan=199; %number of channels
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Loop over all subjects
        cmat=[];
        
        %%%%%%%%%%%%%%%%%%%%%% Create and insert your own data structure input here, uncomment the line below only if you want to see an example of data structure format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %load([path_base '\Inputs\data_example_ft.mat']); % refer to code_for_inputs.m for more details about computation
        load([path_base '\Data\' ftrialname(1).name])
        load([path_base '\Data\' ftrialname(1).name '_source_filter.mat'])
        
        nb_trials = size(data.trial, 2)-1; % because first trial is somehow empty
        data.trial(1) = []; % remove first col
        data.time(1) = []; % remove first col
        
        % Dynamic FC Computation
        cfg = [];
        cfg.window.size = dFC.window.size; % sliding window length in seconds
        cfg.window.step = dFC.window.step; % step between windows in seconds
        cfg.bpfreq = dFC.bpfreq; % frequency band of interest
        cfg.prestim = pre_samples/fs; % prestimulus time in seconds
        cfg.poststim = post_samples/fs; % poststimulus time in seconds
        cfg.conn_method = dFC.conn_method; % FC method
        cfg.labels = scout_labels;
        cmat = go_dynamicFC(cfg,filters,data);
        
        
        save(['F:\WJD\Simon Dynamic FC\Results\FCmat\CAT\' frequency '\winsize_step\' ftrialname(sub_ind).name num2str(dFC.window.size) '_' num2str(dFC.window.step) '_incong.mat'], 'cmat', '-v7.3')
    end
end



%% Calculate number of optimal ICs to search for using DIFFIT

% 1- Define some variables for ICA (To manipulate by the user)

cd([path_base '\Results\FCmat\CAT\' frequency '\winsize_step'])

matname = dir('*incong.mat');


conn_method    = dFC.conn_method; %plv_dyn or wPLI or plv_inst_pn
nb_sub         = length(matname); %number of subjects in the analysis
if strcmp(atlas, 'destrieux')
    nROIs = 148; % 148 destrieux, 68 desikan
elseif strcmp(atlas, 'desikan')
    nROIS = 68;
else
    error('Error in atlas definition')
end

% NCs            = 6; %number of output components (this is a parameter to manipulate with, till now not found a specific automatic method to determine optimal number of independent components)
ICA_val        = 1; % 1=ICA-JADE, 2=ICA-InfoMax, 3=ICA-SOBI, 4=ICA-FastICA, 5=ICA-CoM2, 6:ICA-PSAUD
%%% cfg.val is a parameter to manipulate the subtype of ICA to use in the analysis,
%%% based on our previous analysis, it is preferred to apply ICA that uses high statistical order like JADE,InfoMax,CoM2 and PSAUD
%%% so for exp, we can start by testing ICA-JADE (cfg.val=1)


for i=1:nb_sub
    cmat_list{1} = [matname(i).name];
    %% 2- Load and Concatenate cmat of all subjs and trials for once (to save time and not to repeat for every number of compo)
    cfg.data        = cmat_list;
    cfg.n_parcels   = nROIs;
    [cmat_all,index,sub_trials,cmat_time]=load_and_concatenate(cfg); %all
    %% 3- Calculate DIFFIT values based on 'goodness_of_fit' algorthm to select the optimal number of components for a range of input nb_compo values
    % References:
    % -Tewarie, P., Liuzzi, L., O’Neill, G.C., Quinn, A.J., Griffa, A., Woolrich, M.W., Stam, C.J., Hillebrand, A., Brookes, M.J., 2019b. Tracking dynamic brain networks using high temporal resolution MEG measures of functional connectivity. NeuroImage
    % -Zhu, Y., Liu, J., Ye, C., Mathiak, K., Astikainen, P., Ristaniemi, T., Cong, F., 2020. Discovering dynamic task-modulated functional networks with specific spectral modes using MEG. NeuroImage
    % -Timmerman, M.E., Kiers, H.A., 2000. Three-mode principal components analysis: choosing the numbers of components and sensitivity to local optima. Br J Math Stat Psychol
    % -Wang, D., Zhu, Y., Ristaniemi, T., Cong, F., 2018. Extracting multi-mode ERP features using fifth-order nonnegative tensor decomposition. J Neurosci Methods
    nb_min=3; %the minimum nb of components to evaluate DIFFIT on (minimum allowed value for nb_min is 2)
    nb_max=10; %the maximum nb of components to evaluate DIFFIT on (here it is 12 just an example, however this is to be estimated by the user)
    for nb=nb_min-1:nb_max+1
        cfg = [];
        cfg.NCs =  nb;
        cfg.data = cmat_list;
        cfg.n_parcels = nROIs;
        cfg.val = ICA_val;
        cfg.cmat_all = cmat_all;
        cfg.index = index;
        cfg.sub_trials = sub_trials;
        cfg.cmat_time = cmat_time;
        tic
        results = go_decomposeConnectome_SS_ICA_noload(cfg); %all
        timeElapsed = toc
        mat=results.cmat_all;
        matprime=results.mixing*resultms.sig;
        F(nb,:)=goodnessOfFit_edit(matprime',mat','NRMSE');
        FIT(nb,1)=mean(F(nb,:),2);
    end
    for J=nb_min:nb_max
        DIFFIT(J)=(FIT(J)-FIT(J-1))/(FIT(J+1)-FIT(J));
    end
    [max_DIFFIT,nb_opt]=max(DIFFIT); %nb_opt is considered as the optimal nb of compo based on goodness of fit approach
    GoF(i).F = F;
    GoF(i).FIT = FIT;
    GoF(i).nb_opt = nb_opt;
    GoF(i).max_DIFFIT = max_DIFFIT;
    disp(['subject n ' num2str(i) ' done'])
    clear 'cmat_all' 'index' 'sub_trials' 'cmat_time' 'results' 'mat' 'matprime' 'cfg' 'cmat_2d' 'cmat_red' 'cmat'
end

save(['F:\WJD\Simon Dynamic FC\Results\ICA\HC_PD_CAT\' frequency '\winsize_step\' 'GoF.mat'], 'GoF', '-v7.3')



% Configuration for states (ICA algo)
%% 2- Load and Concatenate cmat of all subjs and trials for once (to save time and not to repeat for every number of compo)
ICA_val        = 1; % 1=ICA-JADE, 2=ICA-InfoMax, 3=ICA-SOBI, 4=ICA-FastICA, 5=ICA-CoM2, 6:ICA-PSAUD

for i=2:length(matname)
    cmat_list{1} = [matname(i).name];
    
    
    cfg.data        = cmat_list;
    cfg.n_parcels   = nROIs;
    
    [cmat_all,index,sub_trials,cmat_time]=load_and_concatenate(cfg); % all
    
    
    %% Run ICA
    cfg = [];
    cfg.NCs = 10;
    cfg.data = cmat_list;
    cfg.n_parcels = nROIs;
    cfg.val = ICA_val;
    cfg.cmat_all = cmat_all;
    cfg.index = index;
    cfg.sub_trials = sub_trials;
    cfg.cmat_time = cmat_time;
    results = go_decomposeConnectome_SS_ICA_noload(cfg); %all
    
    % Permutation testing
    % Configuration for null distribution (Sign-Flip approach)
    
    cfg_null                   = [];
    cfg_null.p                 = 0.05;
    cfg_null.bonferroni_factor = 2*cfg.NCs; % 2 tails x NCs ICs
    cfg_null.nperms =10000; % If nperms not define, default is to take the binomial coeff of (nsub, half nsub) which is very high !
    %perms = go_testNetworks_general(cfg_null, results);
    
    save(['F:\WJD\Simon Dynamic FC\Results\ICA\HC_PD_CAT\' frequency '\' 'winsize_step\' matname(i).name '.mat'], 'results', '-v7.3');
    %save(['F:\WJD\Simon Dynamic FC\Results\ICA\HC_PD_CAT\' frequency '\' 'perms.mat'],'perms');
end
% Extract FC matrices for network plotting

cd('F:\WJD\Simon Dynamic FC\Results\ICA\HC_PD_CAT\beta\winsize_step')
ftrialname = dir('*.mat');
bigmat = zeros(148, 148, 10, length(ftrialname));



corrmat = zeros(3,3,10);
for neti = 1:10
    for rowi = 1:9
        for coli = 1:9
            
            A = reshape(mean(bigmat(:,:,neti,rowi),3), 1, []);
            B = reshape(mean(bigmat(:,:,neti,coli),3), 1, []);
            corrmat(rowi, coli, neti) = corr(A',B');
        end
    end
end

imagesc(squeeze(mean(corrmat, 3)))