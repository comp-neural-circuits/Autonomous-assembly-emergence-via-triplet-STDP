%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used to calculate the graph measures of the manuscript:
% Lisandro Montangie, Christoph Miehl, Julijana Gjorgjieva (2020) PLoS
% Computational Biology
%
%
%
% The following code is used to generate the averaged connectivity matrix
% based on a k-means clustering algorithm.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data_triplet_STDP.mat'); % load data from main simulation

W_clustered=cell(1,length(ETAM)); %Initialize averaged cluster
for gg=1:length(ETAM)
    W_clustered{gg}=zeros(N);
end


% to calculate proportion of bidirectional connections
sum_uni=zeros(rep,length(ETAM));
sum_bid=zeros(rep,length(ETAM));

for jj1=1:length(ETAM)
    for jj2=1:rep
            
        W=Wt{jj1,jj2};
        W=(W>w_max*0.05).*W;

        %% k-means clustering
        k=8; %number of clusters

        [a,b,dis]=kmeans([W W'],k,'EmptyAction','singleton','Replicates',2000);
        % create a permutation matrix that adjust neurons with the same index to
        % the same group:
        [Y,I] = sort(a);
        permu=zeros(N);
        for jj=1:N
            permu(jj,I(jj))=1;
        end
        new_W=permu*W*permu^-1;

        mat1=W(randperm(size(W,1)),:);
        mat_re=mat1(:,randperm(size(W,2)));
        [a2,b2,dis2]=kmeans([mat_re mat_re'],k,'EmptyAction','singleton','Replicates',100);

        % reorder the groups for visualizing the feedforward connectivity:
        % for each group, find to each groups its projects:
        Y_re=zeros(N,1); %reordered  groups indexes vector.
        I_re=zeros(N,1); %reordered the indexes of the neurons.
        % ajust the first group into the first indexes
        n1=numel(find(Y==1));
        Y_re(1:n1)=1*ones(1,n1);
        I_re(1:n1)=I(1:n1);
        used_group=1:k;% the group that were not adjusted yet.
        used_group(1)=0;
        new_ind=1;
        for jj=1:k-1
            % for each group, find for which of the other groups the number of entries is the largest.
            score=zeros(k,1);
            for l=1:k
                % calculate the number of projections into the group l:
                mat1=zeros(N,numel(find(Y==new_ind)));
                mat1(find(Y==l),:)=1;
                score(l)=sum(sum(mat1.*new_W(:,find(Y==new_ind))));
            end
            l_max=find(score==max(score));% the next group
            % check that the next group does not appear already
            if used_group(l_max(1))==0 % if the groups appears already
                %adjust one group that where not used yet
                vec=find(used_group~=0);
                new_ind=vec(1);
            else
                new_ind=l_max(1);
                % adjust group l_max next:

            end
            first_ind=numel(find(Y_re~=0))+1; % the first index of the new group;
            n_l=numel(find(Y==new_ind)); % number of neurons in l_max;
            Y_re(first_ind:first_ind+n_l-1)=new_ind*ones(n_l,1);
            I_re(first_ind:first_ind+n_l-1)=I(Y==new_ind);
            used_group(new_ind)=0;
        end

        % reorder W again:
        permu2=zeros(N);
        for jj=1:N
            permu2(jj,I_re(jj))=1;
        end
        re_W=permu2*W*permu2^-1;        

        W_clustered{1,jj1}=W_clustered{1,jj1}+re_W;
        
        W_tmp=re_W;
        W_tmp=(W_tmp>w_max*0.5);

        % Calculate proportion of bidirectional connections
        sum_total=sum(sum(W_tmp));
        Wun=W_tmp-W_tmp';
        Wbi=W_tmp+W_tmp';  
        Wbid=(Wbi==2);
        Wunid=(Wun==1);  
        sum_uni(jj2,jj1)=sum(sum(Wunid))/sum_total;
        sum_bid(jj2,jj1)=sum(sum(Wbid))/sum_total;
    end
end

for gg=1:length(ETAM)
    W_clustered{gg}=W_clustered{gg}./rep;
    mean_uni(gg)=mean(sum_uni);
    mean_bid(gg)=mean(sum_bid);
    err_uni(gg)=std(sum_uni);
    err_bid(gg)=std(sum_bid);
end


for kk1=1:length(ETAM)
    h3=figure;
    imagesc(W_clustered{kk1})
    colormap(flipud(gray))
    caxis([0,w_max])
    colorbar
    XLABEL=xlabel('# Presynaptic Neuron');
    YLABEL=ylabel('# Postsynaptic Neuron');
    set([XLABEL,YLABEL], ...
        'FontName'   , 'Helvetica');
    set([XLABEL,XLABEL], ...
        'FontSize'   , 16          );
    set(gca,'fontsize',16)
    set(gcf, 'PaperPositionMode', 'auto');
    axis square
    set(h3,'Units','Inches');
    set(gca,'FontSize',16)
    pos = get(h3,'Position');
    set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    name3 = sprintf(['Average_conn_matrix_' num2str(ETAM(kk1)) '.pdf']);
    print(h3,name3,'-dpdf','-r500');   
   
    h4=figure;
    c = categorical({'UC','BC'});
    barwitherr([err_uni(kk1),err_bid(kk1)],[mean_uni(kk1),mean_bid(kk1)],0.4)
    ylim([0 1.09])
    YLABEL=ylabel('Fraction');
    set([XLABEL,YLABEL],'FontName', 'Helvetica');
    set([XLABEL,XLABEL],'FontSize',16);
    set(gca,'fontsize',16)
    set(gcf, 'PaperPositionMode', 'auto');
    set(gca, 'box', 'off')
    set(h4,'Units','Inches');
    set(gca,'FontSize',16)
    pos = get(h4,'Position');
    set(h4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    name4 = sprintf(['barunibid_' num2str(ETAM(kk1)) '.pdf']);
    print(h4,name4,'-dpdf','-r500');
end
