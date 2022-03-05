function go_viewNetworkComponents_eeg_concat_interface(opt,results,varargin)

% This code takes as input opt structure that should contains the basic info for plotting components:
% The results components will be plotted together on a figure using
% subplots (max 6 components per figure), each component is splitted into
% brain network plot (results.maps) and temporal signals plot
% (results.signals averaged over all trials and subjects)

% set some defaults
opt.components = ft_getopt(opt,'components',1:results.NCs);
opt.threshold  = ft_getopt(opt,'threshold',0.9);

if opt.threshold > 1 || opt.threshold < 0
    error('threhsold can only be between 0 and 1');
end

% there may be some statsitical thresholds to overlay, check
if nargin == 3;
    stat_thresh = varargin{1}.thresholds;
end

scout_labels=opt.scout_labels;
scout_mni=opt.scout_mni;
Surfmatrix=opt.Surfmatrix;
meth=opt.meth;

if((floor(length(opt.components)/6)==0)||(length(opt.components)==6))
    figure()
    for ii = opt.components(1):opt.components(length(opt.components))
    mode = results.maps(:,:,ii);
    modescale = max(abs(mode(:)));
    threshmode = mode.*(abs(mode) >= opt.threshold*modescale);
%     if mean(threshmode(:)) < 0;
%         flip = -1;
%     else
%         flip = 1;
%     end
    set(gcf,'Units','normalized')
    set(gcf,'Position',[0.1844    0.4722    0.6401    0.3889]);
    subplot(3,4,2*ii-1)
    go_view_brainnetviewer_eeg_interface(mode,opt.threshold,opt.label,meth,scout_labels,scout_mni,Surfmatrix);
    subplot(3,4,2*ii)
    IC = squeeze(results.signals(ii,:,:));
    if size(IC,1)==1
        plot(results.time,IC,'linewidth',2,'color','k');
    else
        plot(results.time,mean(IC,2),'linewidth',2,'color','k');
    end
    if exist('stat_thresh','var');
        hold on
        fill([results.time fliplr(results.time)],[stat_thresh.upper(ii,:) fliplr(stat_thresh.lower(ii,:))]',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none');
        hold off
    end
     
    grid off
    set(gcf,'color','w')
    set(gca,'xlim',[min(results.time) max(results.time)])
    xlabel('time / s')
    ylabel('trial averaged IC')
 
    end
else
    nb_rep=floor(length(opt.components)/6);
    for nb=1:nb_rep   
        figure()
        for ii = opt.components(1+6*(nb-1)):opt.components(6+6*(nb-1))
            i=ii-6*(nb-1);
            mode = results.maps(:,:,ii);
            modescale = max(abs(mode(:)));
            threshmode = mode.*(abs(mode) >= opt.threshold*modescale);
%             if mean(threshmode(:)) < 0;
%                 flip = -1;
%             else
%                 flip = 1;
%             end
            set(gcf,'Units','normalized')
            set(gcf,'Position',[0.1844    0.4722    0.6401    0.3889]);
            subplot(3,4,2*i-1)
            go_view_brainnetviewer_eeg_interface(mode,opt.threshold,opt.label,meth,scout_labels,scout_mni,Surfmatrix);
            subplot(3,4,2*i)
            IC = squeeze(results.signals(ii,:,:));
            if size(IC,1)==1
                plot(results.time,IC,'linewidth',2,'color','k');
            else
                plot(results.time,mean(IC,2),'linewidth',2,'color','k');
            end
    
            if exist('stat_thresh','var');
                hold on
                fill([results.time fliplr(results.time)],[stat_thresh.upper(ii,:) fliplr(stat_thresh.lower(ii,:))]',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none');
                hold off
            end
            
            grid off
            set(gcf,'color','w')
            set(gca,'xlim',[min(results.time) max(results.time)])
            xlabel('time / s')
            ylabel('trial averaged IC')
        end
    end       
        figure()
        for ii = opt.components(6*nb_rep+1):opt.components(length(opt.components))
            i=ii-6*nb_rep;
            mode = results.maps(:,:,ii);
            modescale = max(abs(mode(:)));
            threshmode = mode.*(abs(mode) >= opt.threshold*modescale);
%             if mean(threshmode(:)) < 0;
%                 flip = -1;
%             else
%                 flip = 1;
%             end
            set(gcf,'Units','normalized')
            set(gcf,'Position',[0.1844    0.4722    0.6401    0.3889]); 
            subplot(3,4,2*i-1)
            go_view_brainnetviewer_eeg_interface(mode,opt.threshold,opt.label,meth,scout_labels,scout_mni,Surfmatrix);
            subplot(3,4,2*i)
            IC = squeeze(results.signals(ii,:,:));
            if size(IC,1)==1
                plot(results.time,IC,'linewidth',2,'color','k');
            else
                plot(results.time,mean(IC,2),'linewidth',2,'color','k');
            end
    
            if exist('stat_thresh','var');
                hold on
                fill([results.time fliplr(results.time)],[stat_thresh.upper(ii,:) fliplr(stat_thresh.lower(ii,:))]',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none');
                hold off
            end
     
            grid off
            set(gcf,'color','w')
            set(gca,'xlim',[min(results.time) max(results.time)])
            xlabel('time / s')
            ylabel('trial averaged IC')
        end
end
end