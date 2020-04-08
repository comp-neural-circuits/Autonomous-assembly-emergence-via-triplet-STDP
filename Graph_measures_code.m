%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This code is used to calculate the graph measures of the manuscript:
% Lisandro Montangie, Christoph Miehl, Julijana Gjorgjieva (2020) PLoS
% Computational Biology
%
% The graph measures are from the Brain Connectivity Toolbox 
% (brain-connectivity-toolbox.net): Rubinov M, Sporns O (2010) NeuroImage
% 52:1059-69. 
%
% The following code is used to generate the three graph measurements as in
% Fig. 7-11, Suppl. Fig. 2-4.
%
%
% barwitherr function from:
% Martina Callaghan (2020). barwitherr(errors,varargin) 
% (https://www.mathworks.com/matlabcentral/fileexchange/30639-barwitherr-errors-varargin),
% MATLAB Central File Exchange. Retrieved March 16, 2020. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

load('data_triplet_STDP.mat');

% Add the BCT folder to path
addpath(genpath([pwd,'\BCT']))

xaxis=ETAM;
%xaxis=TAUE;

% Initialize matrices
meancc=zeros(length(xaxis),1);
meaneg=zeros(length(xaxis),1);
meanmod=zeros(length(xaxis),1);
errcc=zeros(length(xaxis),1);
erreg=zeros(length(xaxis),1);
errmod=zeros(length(xaxis),1);

for jj1=1:length(xaxis)

    ccr=0;
    egr=0;
    Qr=0;

    for jj2=1:rep
        W=(Wt{jj1,jj2}>w_max*0.05).*Wt{jj1,jj2};

        cc=clustering_coef_wd(W); % Calculate clustering coefficient
        Eglob = efficiency_wei(W); % Calculate global efficiency
        [Ci,Q]=modularity_dir(W,1); % Calculate modularity

        ccr=[ccr;mean(cc)];
        egr=[egr;Eglob];
        Qr=[Qr;Q];

    end
    ccr(1)=[];
    egr(1)=[];
    Qr(1)=[];

    meancc(jj1)=mean(ccr);
    meaneg(jj1)=mean(egr);
    meanmod(jj1)=mean(Qr);

    errcc(jj1)=3*std(ccr)/sqrt(rep);
    erreg(jj1)=3*std(egr)/sqrt(rep);
    errmod(jj1)=3*std(Qr)/sqrt(rep);

end

S{1}.mcc=meancc;
S{1}.mge=meaneg;
S{1}.mm=meanmod;
S{1}.ecc=errcc;
S{1}.ege=erreg;
S{1}.em=errmod;


h1=figure
errorbar(ETAM,S{1}.mcc,S{1}.ecc,'LineWidth',3);
XLABEL=xlabel('Modulation parameter \eta_-');
YLABEL=ylabel('Clustering coefficient');
set([XLABEL,YLABEL],'FontName','Helvetica');
set([XLABEL,XLABEL],'FontSize',16);
set(gca,'FontSize',16,'FontName','Helvetica');
set(gca,'linewidth',3)
set(h1,'Units','Inches');
pos1 = get(h1,'Position');
set(gca, 'box', 'off');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos1(3), pos1(4)]);
name1 = sprintf('clustcoeff.pdf');
print(h1,name1,'-dpdf','-r500');  

h2=figure
errorbar(ETAM,S{1}.mge,S{1}.ege,'LineWidth',3);
XLABEL=xlabel('Modulation parameter \eta_-');
YLABEL=ylabel('Global efficiency');
set([XLABEL,YLABEL],'FontName','Helvetica');
set([XLABEL,XLABEL],'FontSize',16);
set(gca,'FontSize',16,'FontName','Helvetica');
set(h2,'Units','Inches');
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'linewidth',3)
pos2 = get(h2,'Position');
set(gca, 'box', 'off')
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos2(3), pos2(4)]);
name2 = sprintf('globaleff.pdf');
print(h2,name2,'-dpdf','-r500');   


h3=figure
errorbar(ETAM,S{1}.mm,S{1}.em,'LineWidth',3);hold on;
XLABEL=xlabel('Modulation parameter \eta_-');
YLABEL=ylabel('Modularity');
set([XLABEL,YLABEL],'FontName','Helvetica');
set([XLABEL,XLABEL],'FontSize',16);
set(gca,'fontsize',16)
set(h3,'Units','Inches');
set(gca,'FontSize',16,'FontName','Helvetica');
set(gca,'linewidth',3)
pos3 = get(h3,'Position');
set(gca, 'box', 'off')
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos3(3), pos3(4)]);
name3 = sprintf('modularity.pdf');
print(h3,name3,'-dpdf','-r500');   

