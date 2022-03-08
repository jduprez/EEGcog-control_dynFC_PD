function [thresholded_mat]=threshold_FCmat(C,thresh,meth)

% This function thresholds the FC matrix according to different methods.
% Possible methods :
% - thresh_abs; threshold=0.93 means that we want to display all nodes with connections value higher than 0.93*max(connectivity values) ==> so components will have varied number of nodes depending on connectivities values
% - thresh_val; threshold=0.93 means that we want to display all nodes having a strength higher than 0.93*max(all node_strength values) ==> so components will have varied number of nodes depending on node strength values
% - thresh_node; threshold=0.93 means that we want to display the (floor(0.93*total_nb_nodes)) nodes that have the highest node strength values ==> so all components will have this fixed number of nodes with the strongest node strength values
% - thresh_conn; threshold=0.93 means that we want to display the top
% x=floor((1-0.93)*(nROI*(nROI-1)/2)) edges with the highest FC values

nROI=size(C,1);
%nodes strength
for r=1:nROI
    nodeStrength(r)=sum(abs(C(:,r)));
end

if(strcmp(meth,'thresh_abs'))
    edge_thresh=C;
    C(eye(size(C))==1) = 0;
    limit = max(abs(C(:)))*thresh;
    % C=abs(C);
    mask = abs(C) >= limit;
    cLims = [-max(abs(C(:))) max(abs(C(:)))];
    % cLims = [0 max(abs(C(:)))];
    C(~mask) = NaN;
    edge_thresh(~mask)=0;
elseif(strcmp(meth,'thresh_val'))
    cLims = [-max(abs(C(:))) max(abs(C(:)))];
    str_thresh_val=thresh*max(nodeStrength);
    str_thresh_ind=find(nodeStrength>str_thresh_val);
    C_thresh=zeros(nROI,nROI);
    for j=1:nROI
        for k=1:nROI
            if(ismember(j,str_thresh_ind)&&ismember(k,str_thresh_ind))
                C_thresh(j,k)=C(j,k);
            end
        end
    end
    C_thresh(C_thresh==0)=NaN;
    C=C_thresh;
elseif(strcmp(meth,'thresh_node'))
    cLims = [-max(abs(C(:))) max(abs(C(:)))];
%     [~,I] = maxk(nodeStrength,1+floor((1-thresh)*nROI));
    [~,ind]=sort(nodeStrength,'descend');
    I=ind(1:1+floor((1-thresh)*nROI));
    C_thresh=zeros(nROI,nROI);
    for j=1:nROI
        for k=1:nROI
            if(ismember(j,I)&&ismember(k,I))
                C_thresh(j,k)=C(j,k);
            end
        end
    end
    C_thresh(C_thresh==0)=NaN;
    C=C_thresh;
    
    elseif(strcmp(meth,'thresh_conn'))
    cLims = [-max(abs(C(:))) max(abs(C(:)))];
    Cu=triu(C);
    max_Cu=maxk(abs(Cu(:)),floor((1-thresh)*(nROI*(nROI-1)/2)));
    C_thresh=zeros(nROI,nROI);
    for kk=1:length(max_Cu)
         if(sum(any(Cu==max_Cu(kk))))
            [rr(kk),cc(kk)]=find((Cu==max_Cu(kk)));
        elseif(sum(any(Cu==-max_Cu(kk))))
            [rr(kk),cc(kk)]=find((Cu==-max_Cu(kk)));
        end        
        C_thresh(rr(kk),cc(kk))=Cu(rr(kk),cc(kk));
    end
    C_thresh=C_thresh+C_thresh';
    C_thresh(C_thresh==0)=NaN;
    thresholded_mat=C_thresh;
end

end