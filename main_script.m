%% Info TBD

%% 1.1. Define basic param and add basic folders to matlab path

path_base = 'F:\WJD\Simon Dynamic FC';
mri_template = 'icbm'; % 'colin' for Colin27 or 'icbm' for ICBM152
atlas = 'destrieux'; % 'desikan' for Desikan68 atlas or 'destrieux' for Destrieux148 atlas

% Automatically add folders and subfolders (Code and Inputs) necessary for
% loading variables and executing codes
addpath(genpath(['F:\WJD\Simon Dynamic FC\input4code\' mri_template]));
addpath(genpath('C:\GitHub\EEGcog-control_dynFC_PD'));

% Add openmeeg path
setenv('PATH', [path_base '\input4code\Toolboxes\OpenMEEG\bin'])
addpath([path_base '\input4code\Toolboxes\fieldtrip-20190224\external\openmeeg'])

% Add default fieldtrip path
tmp = which('ft_defaults');
if isempty(tmp)
    addpath([path_base '\input4code\Toolboxes\fieldtrip-20190224']); %add defaults fieldtrip
    ft_defaults
end

%% 1.2. Define wMNE specification (prefer to keep these default params values)

source.weightExp   = 0.5; % param for depth weighting
source.weightLimit = 10; % param for depth weighting limit
source.SNR         = 3; % signal to Noise Ratio

%% 1.3. Define dynamic Functional Connectivity (dFC) specification

dFC.conn_method = 'plv_dyn'; %'plv_dyn' for windowedPLV or 'wPLI' for windowedwPLI or 'plv_inst_pn' for instantaneous PLV
dFC.bpfreq      = [30 45]; %frequency band of interest
dFC.window.size = 0.16; % sliding window length in seconds (for example calculated for 6cycles,CentralFreq=35 ==> 6/35=0.17s), in case of 'plv_inst_pn' this input is meaningless
dFC.window.step = 0.016; % step between windows in seconds (for example 90% overlapping=10/100*window_size), in case of 'plv_inst_pn' this input is meaningless   
% dFC.bpfreq      = [12 25]; %frequency band of interest
% dFC.window.size = 0.3243; % sliding window length in seconds (for example calculated for 6cycles,CentralFreq=35 ==> 6/35=0.17s), in case of 'plv_inst_pn' this input is meaningless
% dFC.window.step = 0.03243; % step between windows in seconds (for example 90% overlapping=10/100*window_size), in case of 'plv_inst_pn' this input is meaningless   


%% 2.1. Read MRI download from freesurfer (Colin27 or ICBM152) + Realign MRI + Compute Hheadmodel (OpenMEEG)

% Run this section once (as the headmodel is common between all subjects
% when we are using template MRI): Output: subvol (=headmodel). Can be saved and used after loading instead of computing several times.

%load realigned MRI saved in Inputs
load([path_base '\input4code\' mri_template '\mri_' mri_template '_realign.mat']); 

% BEM Headmodel: Compute the BEM headmodel using OpenMEEG

cfg           = [];
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri_realign); %ctf/mm

cfg = [];
cfg.method = 'headshape';
brain = ft_read_headshape([path_base '\Inputs\' mri_template '\tess_innerskull_' mri_template '.mat']); % load the innerskull template used in Brainstorm 
brain_mm = ft_convert_units(brain,'mm');
cfg.headshape = brain_mm;
cfg.numvertices = [3000];
bnd(1) = ft_prepare_mesh(cfg,segmentedmri);

cfg = [];
cfg.method = 'headshape';
skull = ft_read_headshape([path_base '\Inputs\' mri_template '\tess_outerskull_' mri_template '.mat']); % load the outerskull template used in Brainstorm (colin27)
skull_mm = ft_convert_units(skull,'mm');
cfg.headshape = skull_mm;
cfg.numvertices = [3000];
bnd(2) = ft_prepare_mesh(cfg,segmentedmri);

cfg = [];
cfg.method = 'headshape';
head = ft_read_headshape([path_base '\Inputs\' mri_template '\tess_head_' mri_template '.mat']); % load the head template used in Brainstorm (colin27)
head_mm = ft_convert_units(head,'mm');
cfg.headshape = head_mm;
cfg.numvertices = [3000];
bnd(3) = ft_prepare_mesh(cfg,segmentedmri);

figure();
ft_plot_mesh(bnd(1), 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(bnd(2), 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(bnd(3), 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3

cfg = [];
cfg.method ='openmeeg'; 
subvol = ft_prepare_headmodel(cfg, bnd); % ctf/mm

%% 2.2 Load Some Variables for further usage (Run the section as it is)

% Load already computed and saved variables (as we are dealing with template
% MRI; these variables were computed for once as in code_for_inputs.m function
% and loaded here to be used for all subjects

load([path_base '\input4code\' mri_template '\elec_BS_' mri_template '_199ch.mat']); % fieldtrip format electrode already computed from modified EGI257 system: extracted from BrainStorm
load([path_base '\input4code\' mri_template '\' atlas '\' atlas '_scs_' mri_template '.mat']); % scout structure to extract sources centroids positions orientations (in SCS or CTF coord): calculated in Section 2 of 'code_for_inputs'
load([path_base '\input4code\' mri_template '\' atlas '\subgrid.mat']); % fieldtrip format of the grid specific to the mri and scout used (refer to code_for_inputs.m for details about computation): calculated in Section 3 of 'code_for_inputs'
load([path_base '\input4code\' atlas '_labels.mat']); % cell labels for the scout regions: extracted from BrainStorm

%% 3. Main code to compute sources+connectivities

cd('D:\JoanD\EEGCOG_data_4_DynEEG')
ftrialname = dir('*_cong.mat');

%%%%%%%%%%%%%%%%%%%%%% Manipulate these parameters appropriatly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_sub=length(ftrialname); %number of subjects
%nb_trials=200; %number of trials per subject (if variable between subjects, define this param inside the for loop subjects)
fs=1000; %sampling frequency
pre_samples=700; %nb samples before trial onset 
post_samples=1200; %nb samples after trial onset
nb_samples=1901; %total number of samples in each trial (again, if variable between subjects, define this param inside the for loop subjects, preferable to be consistent between trials for a single subject for later concatenation)
nb_chan=199; %number of channels (electrodes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Loop over all subjects
for sub_ind = 20:nb_sub
%     filters=[];
    cmat=[];
    
    %%%%%%%%%%%%%%%%%%%%%% Create and insert your own data structure input here, uncomment the line below only if you want to see an example of data structure format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %load([path_base '\Inputs\data_example_ft.mat']); % refer to code_for_inputs.m for more details about computation
    load(['D:\JoanD\EEGCOG_data_4_DynEEG\' ftrialname(sub_ind).name])
    %load(['D:\JoanD\EEGCOG_data_4_DynEEG\' ftrialname(sub_ind).name '_source_filter.mat'])
    nb_trials = size(data.trial, 2)-1; % because first trial is somehow empty
    data.trial(1) = []; % remove first col
    data.time(1) = []; % remove first col
    % Create your own data structure with the following fields: 
    % data.label: cell of length nb_chan, each cell is a string of channels or electrodes labels (for example in case of EGI257, we can label electrodes from E1 to E256 and Cz for the last electrodes as a reference;
    % data.fsample: int indicating sampling frequency of data, in our case =fs;
    % data.trial: the segmented data, cell of length nb_trials, each cell is 2D matrix [nb_chan*nb_samples];
    % data.time: the time in sec from -prestim to +poststim, cell of length nb_trials, each cell is 1D matrix [1*nb_samples];
    % data.elec: fieldtrip structure of electrodes, in our case =elec_BS_mm;
    
    % NOTE: DATA TRIALS SHOULD BE PREPROCESSED, SEGMENTEND RELATIVE TO ONSET TRIAL AND CENTRALIZED TO ZERO
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Source Reconstruction
    cfg             = [];
    cfg.prestim     = pre_samples/fs; %prestimulus in seconds for noise cov computation
    cfg.weightExp   = source.weightExp; %param for depth weighting
    cfg.weightLimit = source.weightLimit; %param for depth weighting limit
    cfg.SNR         = source.SNR; %Signal to Noise Ratio
    filters         = go_source_reconstruction(cfg,data,subvol,subgrid,scout_scs.orients);

    %Dynamic FC Computation
    cfg             = [];
    cfg.window.size = dFC.window.size; %sliding window length in seconds 
    cfg.window.step = dFC.window.step; %step between windows in seconds 
    cfg.bpfreq      = dFC.bpfreq; %frequency band of interest
    cfg.prestim     = pre_samples/fs; %prestimulus time in seconds
    cfg.poststim    = post_samples/fs; %poststimulus time in seconds
    cfg.conn_method = dFC.conn_method; %3 choices: 'plv_inst_pn' for instantaneousPLV or 'plv_dyn' for windowedPLV or 'wPLI' for windowedwPLI
    cfg.labels      = scout_labels;
    cmat            = go_dynamicFC(cfg,filters,data);
    
%    save([path_base '\Results\' cfg.conn_method '_cmat_' int2str(sub_ind) '.mat'],'cmat');

    %Uncomment below only if you want to plot connectivity network results
%     %at specific time interval
%     load([path_base '\Inputs\' mri_template '\Surfmatrix_' mri_template '.mat']); %structure extracted from brainnetviewer to draw the brain cortex
%     load([path_base '\Inputs\' mri_template '\' atlas '\' atlas '_mni_' mri_template '.mat']); %scout structure to extract sources centroids positions orientations (in MNI coord)
%     load([path_base '\Inputs\' atlas '_labels.mat']); %cell labels for the scout regions
%     BSN1=cmat.connectivity{1}(:,:,20); %example connectivity matrix at window 20 in case of windowed versions or sample 20 in case pf PLVinstantaneous
% %     BSN1=squeeze(mean(cmat.connectivity{1}(:,:,307:409),3));
%     meth='thresh_node'; %'thresh_abs' for thresholding on absolute values, 'thresh_val' for thresholding on node strength values, 'thresh_node' for thresholding on node strength density
%     thresh_val=0.85; %if 0.85 on thresh_node: meaning that we want to retain the 15% of total number of nodes that have the highest node strength value
%     label_id=0; %0 if we want to plot networks without ROIs labels text, 1 with ROIs labels text
%     figure();
%     go_view_brainnetviewer_eeg_interface(BSN1,thresh_val,label_id,meth,scout_labels,scout_mni,Surfmatrix);
save(['D:\JoanD\EEGCOG_data_4_DynEEG_FCmat\' ftrialname(sub_ind).name '_source_filter.mat'], 'filters', '-v7.3')

save(['D:\JoanD\EEGCOG_data_4_DynEEG_FCmat\' ftrialname(sub_ind).name '_cong_wplv3045.mat'], 'cmat', '-v7.3')

disp(['subject ' num2str(sub_ind) ' is done'])
end

