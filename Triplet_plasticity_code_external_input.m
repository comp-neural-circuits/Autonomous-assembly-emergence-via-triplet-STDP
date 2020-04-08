%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is the basis of the manuscript:
% Lisandro Montangie, Christoph Miehl, Julijana Gjorgjieva (2020) PLoS
% Computational Biology
%
%
% The following code is the basis to investigate the emergence of
% assemblies in the presence of structured external input.
% The code is based on code from Ravid Tannenbaum, Burak (2016) PLoS Comput Biol 12(8):e1005056
%
%
% With this code, Fig.10 from the manuscript can be
% generated. For a detailed explanation see the README file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STDP parameters
Am=0.01; % LTD window amplitude
taup=16.8; % LTP time constant (ms)
tauy=114; % Second LTP time constant (ms)
taum=33.7; % LTD time constant (ms)
taue=5; % First membrane time constant (ms)
taue2=5; % Second membrane time constant (ms)
tauef=taue.*taue2./(taue+taue2);
tauef2=taue.*taue2./(taue+2.*taue2);
eps_0=1; % Normalization of EPSP kernel

%% Simulation parameters
etam=1; % Depression modulation parameter
etap=1; % Potentiation modulation parameter
etay=1; % tauy moduluation parameter

N=48; % Number of neurons in the network
rep=100; % Number of repetitions
mu0=.15; % Input firing rate (ms^(-1))
b_0=ones(N,1)*mu0;
iterations=300000; % Number of time-steps

% Heterosynaptic competition parameters
M=5; %Number of strong synapses
w_max=0.85/M; % keeps the maximum sum of each row and column synaptic efficacy constant so that spectral radius <1 and approximation is good
W_max=M*w_max;% maximum row and coloumn synaptic strength
Psi=0.1/(M*(w_max.^2)); % heterosynaptic competition rate
nu=Psi*0.0005; % learning rate


% External correlated input
CORRINP=[0.25,0.75].*w_max*0.5;

Wt=cell(length(CORRINP),rep);

savefile='data_triplet_STDP';

for kk= 1:numel(CORRINP)
    corr=CORRINP(kk);
    extcorr=corr.*blkdiag(ones(6)-eye(6),ones(6)-eye(6),ones(6)-eye(6),ones(6)-eye(6),ones(6)-eye(6),ones(6)-eye(6),ones(6)-eye(6),ones(6)-eye(6));
    
    % Calculation of motif coefficients (for details see Methods and Supplementary information 2)
    % Figure 5 shows all motifs used in this simulation (except the zero-order motif if delta is non-zero)
    
    % motif alpha=1, beta=0
    M10=2*pi.*Am*taum.*((etap.*taup)./((taue+etap.*taup).*(tauef+etap.*taup))+(etay.*tauy).^2./((taue+etay.*tauy).*(etap.*taup+etay.*tauy).*(tauef+etay.*tauy)));
    
    % motif alpha=0, beta=1
    M01=-2*pi.*Am*taum.*((etam.*taum)./((taue+etam.*taum).*(tauef+etam.*taum))-(etap.*taup).^2./((taue+etap.*taup).*(etap.*taup+etay.*tauy).*(tauef+etap.*taup)));
    
    % motif alpha=1, beta=1
    M11=pi.*Am*taum.*((etap.*taup+tauef2)./((taue+etap.*taup).*(tauef+etap.*taup))+2.*(etap.*taup).^4./((etap.*taup+etay.*tauy).*((etap.*taup).^2-taue.^2).*((etap.*taup).^2-tauef.^2)) ...
        +(taue+taue2).^2./((taue+2.*taue).*(taue-etap.*taup).*(taue+etay.*tauy))+taue2.^3./((taue+taue2).*(etap.*taup-tauef).*(etay.*tauy+tauef)) ...
        -(etam.*taum+tauef2)./((taue+etam.*taum).*(tauef+etam.*taum)));
    
    % motif alpha=2, beta=0
    M20=2*pi.*Am*taum.*((etap.*taup).^3./((taue+etap.*taup).^2.*(tauef+etap.*taup).^2)+(etay.*tauy).^4./((taue+etay.*tauy).^2.*(etap.*taup+etay.*tauy).*(tauef+etay.*tauy).^2));
    
    % motif alpha=0, beta=2
    M02=-2*pi.*Am*taum.*((etam.*taum).^3./((taue+etam.*taum).^2.*(tauef+etam.*taum).^2)-(etap.*taup).^4./((taue+etap.*taup).^2.*(etap.*taup+etay.*tauy).*(tauef+etap.*taup).^2));
    
    % motif alpha=2, beta=1
    M21=pi/2.*Am.*taum.*((13*(etap.*taup.*taue).^3+taue.^4.*(taue2+etap.*taup).^2+taue.^3.*(taue2+etap.*taup).*(2*taue2+etap.*taup).*(taue2+3.*etap.*taup)...
        +etap.*taup.*taue.^2.*taue2.*(3.*taue2+4.*etap.*taup).^2+2.*(etap.*taup.*taue2).^2.*taue.*(8*taue2+13*etap.*taup))./((taue+taue2).*(taue+2.*taue2).^2.*(taue+etap.*taup).^2.*(tauef+etap.*taup).^2)...
        + (taue + taue2).*(etay.*tauy.*taue.^2.*taue2.*(taue.^3.*taue2.^2.*(taue+2.*taue2)+etay.*tauy.*taue.^2.*taue2.*(2*taue.^2+9*taue.*taue2+9*taue2.^2)+(etay.*tauy).^3.*(3*taue.^3+16*taue.^2.*taue2+26*taue.*taue2.^2 ...
        +13*taue2.^3)+(etay.*tauy).^2.*taue.*(taue.^3+10*taue.^2.*taue2+24.*taue.*taue2.^2+16.*taue2.^3))+etap.*taup.*taue.*(taue+2.*taue2).*(taue.^4.*taue2.^3+3.*etay.*tauy.*taue.^3.*taue2.^2.*(taue+2.*taue2)...
        +3.*(etay.*tauy).^2.*taue.^2.*taue2.*(taue.^2+5*taue.*taue2+5*taue2.^2)+(etay.*tauy).^4.*(3.*taue.^3+17*taue.^2.*taue2+28.*taue.*taue2.^2+14*taue2.^3)+(etay.*tauy).^3.*taue.*(taue.^3+12.*taue.^2.*taue2...
        +30*taue.*taue2.^2+20*taue2.^3))+(etap.*taup).^2.*(4*(etay.*tauy).^4.*(taue.^2+3.*taue.*taue2+2*taue2.^2).^2+taue.^4.*taue2.^2.*(taue.^2+3*taue.*taue2+3*taue2.^2)+etay.*tauy.*taue.^3.*taue2.*(2*taue.^3+11*taue.^2.*taue2...
        +21*taue.*taue2.^2 +14*taue2.^3)+(etay.*tauy).^2.*taue.^2.*(taue.^4+11*taue.^3.*taue2+38*taue.^2.*taue2.^2+54.*taue.*taue2.^3+27.*taue2.^4)+(etay.*tauy).^3.*taue.*(3*taue.^4+23*taue.^3.*taue2...
        +62.*taue.^2.*taue2.^2+70.*taue.*taue2.^3+28.*taue2.^4)))./((etap.*taup + etay.*tauy).*(etap.*taup+taue).*(etay.*tauy+taue).^2.*(taue+2.*taue2).^2.*(taue.*taue2+etap.*taup.*(taue+taue2)).*(taue.*taue2+etay.*tauy.*(taue+taue2)).^2)...
        -((taue+taue2).^3./(taue.*(taue+2.*taue2).^3.*(taue+etam.*taum))-taue2.^3./(taue.*(taue+2.*taue2).^2.*(tauef+etam.*taum))));

    % motif alpha=1, beta=2
    M12=-pi/2.*Am.*taum.*((13*(etam.*taum.*taue).^3+taue.^4.*(taue2+etam.*taum).^2+taue.^3.*(taue2+etam.*taum).*(2*taue2+etam.*taum).*(taue2+3.*etam.*taum)...
        +etam.*taum.*taue.^2.*taue2.*(3.*taue2+4.*etam.*taum).^2+2.*(etam.*taum.*taue2).^2.*taue.*(8*taue2+13*etam.*taum))./((taue+taue2).*(taue+2.*taue2).^2.*(taue+etam.*taum).^2.*(tauef+etam.*taum).^2)...
        -(taue+taue2).*(etay.*tauy.*taue.^4.*taue2.^2.*(taue.*taue2.*(taue+2.*taue2)+etay.*tauy.*(taue.^2+3*taue.*taue2+3*taue2.^2))+etap.*taup.*taue.^3.*taue2.*(taue+2*taue2).*(taue.^2.*taue2.^2 ...
        +3*etay.*tauy.*taue.*taue2.*(taue+2*taue2)+(etay.*tauy).^2.*(2*taue.^2+7*taue.*taue2+7*taue2.^2))+(etap.*taup).^3.*taue.*(taue+2*taue2).*(taue.^2.*taue2.*(taue.^2+8.*taue.*taue2+8*taue2.^2) ...
        +(etay.*tauy).^2.*(3*taue.^3+17*taue.^2.*taue2+28*taue.*taue2.^2+14*taue2.^3)+etay.*tauy.*taue.*(taue.^3+12.*taue.^2.*taue2+30*taue.*taue2.^2+20*taue2.^3)) ...
        +(etap.*taup).^4.*(taue+taue2).*(4*(etay.*tauy).^2.*(taue+taue2).*(taue+2*taue2).^2+taue.^2.*taue2.*(3*taue.^2+13*taue.*taue2+13*taue2.^2)+etay.*tauy.*taue.*(3*taue.^3+20.*taue.^2.*taue2 ...
        +42*taue.*taue2.^2+28*taue2.^3))+(etap.*taup).^2.*taue.^2.*(taue.^2.*taue2.^2.*(2*taue.^2+9*taue.*taue2+9*taue2.^2)+3*etay.*tauy.*taue.*taue2.*(taue.^3+7*taue.^2.*taue2+15*taue.*taue2.^2+10*taue2.^3) ...
        +(etay.*tauy).^2.*(taue.^4+11*taue.^3.*taue2+38*taue.^2.*taue2.^2+54*taue.*taue2.^3+27*taue2.^4)))./((etap.*taup+etay.*tauy).*(etap.*taup+taue).^2.*(etay.*tauy+taue).*(taue+2*taue2).^2.*(taue.*taue2 ...
        +etap.*taup.*(taue+taue2)).^2.*(taue.*taue2+etay.*tauy.*(taue+taue2)))...
        -((taue+taue2).^3./(taue.*(taue+2.*taue2).^3.*(taue+etap.*taup))-taue2.^3./(taue.*(taue+2.*taue2).^2.*(tauef+etap.*taup))));
    
    % motif alpha=3, beta=0
    M30=2*pi.*Am*taum.*((etap.*taup).^5./((taue+etap.*taup).^3.*(tauef+etap.*taup).^3)+(etay.*tauy).^6./((taue+etay.*tauy).^3.*(etap.*taup+etay.*tauy).*(tauef+etay.*tauy).^3));
    
    % motif alpha=0, beta=3
    M03=-2*pi.*Am*taum.*((etam.*taum).^5./((taue+etam.*taum).^3.*(tauef+etam.*taum).^3)-(etap.*taup).^6./((taue+etap.*taup).^3.*(etap.*taup+etay.*tauy).*(tauef+etap.*taup).^3));
    
    % motif alpha=2, gamma=0 from auto-covariance C_ii
    Ma2g0=2*pi.*Am*taum.*((etay.*tauy).^3./((taue+etay.*tauy).^2.*(tauef+etay.*tauy).^2)); %the rate is considered in the sum for this terms
    
    % motif alpha=3, gamma=0 from auto-covariance C_ii
    Ma3g0=2*pi.*Am*taum.*((etay.*tauy).^5./((taue+etay.*tauy).^3.*(tauef+etay.*tauy).^3)); %the rate is considered in the sum for this terms
    %the rate is considered in the sum for this terms
    
    % motif alpha=1, gamma=1 from auto-covariance C_ii
    Ma1g1=pi.*Am*taum.*((etay.*tauy+tauef2)./((taue+etay.*tauy).*(tauef+etay.*tauy)));
    
    % motif alpha=2, gamma=1 and alpha=1, gamma=2 from auto-covariance C_ii
    Ma2g1plusa1g2=pi/2.*Am.*taum.*((taue+taue2).*(2*taue.^4.*(taue2+etay.*tauy).^2+etay.^2.*taue2.^3.*tauy.^2.*(-5*etap.*taup+21*etay.*tauy)...
        +2*etay.*taue.*taue2.^2.*tauy.*(7*etay.*tauy.*(2*taue2+3*etay.*tauy)-etap.*taup.*(2*taue2+5.*etay.*tauy))...
        +taue.^3.*(taue2+etay.*tauy).*(4*taue2.^2+etay.*tauy.*(-etap.*taup+5*etay.*tauy)+taue2.*(-etap.*taup+13*etay.*tauy))...
        +taue.^2.*taue2.*(-etap.*taup.*(taue2.^2+6*etay.*taue2.*tauy+6*etay.^2.*tauy.^2)+etay.*tauy.*(17*taue2.^2+42*etay.*taue2.*tauy...
        +26*etay.^2.*tauy.^2))))./((taue+2*taue2).^2.*(taue+etay.*tauy).^2.*(etay.*taue2.*tauy+taue.*(taue2+etay.*tauy)).^2);
    
    % motif alpha=1, beta=0, gamma=1 from third order cumulant K_ij
    Ma1b0g1=4*pi^2.*Am.*taum.* ((etap.*taup).^3.*(taue+taue2).^2.*(3*etay.*tauy.*taue.*taue2+2*etap.*taup.*(taue*taue2+etay.*tauy.*(taue+2*taue2))))...
        ./((2*etap.*taup+taue).*(etay.*tauy.*taue+etap.*taup.*(etay.*tauy+taue)).*(etap.*taup.*etay.*tauy.*taue+etay.*tauy.*taue.*taue2...
        +etap.*taup.*(etay.*tauy+taue).*taue2).*(taue.*taue2+2.*etap.*taup.*(taue+taue2)).*(taue.*taue2+etap.*taup.*(taue+2*taue2)));
    
    % motif alpha=1, beta=1, gamma=1 from third order cumulant K_ij
         Ma1b1g1=(2/3)*pi^2*Am.*taum.*(taue+taue2).^3.*(1/(taue.^2.*(2*etap.*taup+taue).*(-etay.*tauy+taue).*(taue+3.*taue2))+taue2.^3./(taue.^2.*(taue+taue2).*(2*taue+3*taue2).*(taue.*taue2 ...
             +2*etap.*taup.*(taue + taue2)).*(taue.*taue2-etay.*tauy.*(taue+taue2)))-(3.*taue2.^2.*(-2*taue.*taue2+etay.*tauy.*(taue+2*taue2)))./(taue.^2.*(-etay.*tauy+taue).*(taue ...
             +3*taue2).*(2*taue+3*taue2).*(-taue.*taue2+etay.*tauy.*(taue+taue2)).*(taue.*taue2+etap.*taup.*(taue+2*taue2)))-(2*etay.*tauy.*taue+7*etay.*tauy.*taue2+2*taue.*taue2) ... 
             ./((etap.*taup-taue).*taue.*(2*etay.*tauy+taue).*(taue+3*taue2).*(2*taue+3*taue2).*(taue.*taue2+etay.*tauy.*(taue+2*taue2)))+(3.*(etay.*tauy).^4)./((etay.*tauy-taue)...
             .*taue.*(2*etay.*tauy+taue).*(etay.*tauy.*taue+etap.*taup.*(etay.*tauy+taue)).*(-taue.*taue2+etay.*tauy.*(taue+taue2)).*(taue.*taue2+etay.*tauy.*(taue+2*taue2))) ...
             -(taue2.^3.*(5*etay.*tauy.*taue+7.*etay.*tauy.*taue2+2*taue.*taue2))./(taue.*(taue+taue2).*(taue+3*taue2).*(2*taue+3*taue2).*(taue.*taue2-etap.*taup.*(taue+taue2)) ...
             .*(taue.*taue2+2*etay.*tauy.*(taue+taue2)).*(taue.*taue2+etay.*tauy.*(taue+2*taue2)))-(3*(etay.*tauy).^4.*taue2.^2)./((etay.*tauy-taue).*taue.*(-taue.*taue2+etay.*tauy.*(taue+taue2)).*(taue.*taue2 ...
             +2*etay.*tauy*(taue+taue2)).*(taue.*taue2+etay.*tauy.*(taue+2*taue2)).*(etay.*tauy.*taue.*taue2+etap.*taup.*(taue.*taue2+etay.*tauy.*(taue+taue2))))+(3*(etap.*taup).^5.*(3*etay.*tauy.*taue.*taue2 ...
             +2*etap.*taup*(taue.*taue2+etay.*tauy.*(taue+2*taue2))))./((etap.*taup-taue).*(2*etap.*taup+taue).*(etay.*tauy.*taue+etap.*taup.*(etay.*tauy+taue)).*(-taue.*taue2+etap.*taup.*(taue+taue2))...
             .*(taue.*taue2+2*etap.*taup.*(taue+taue2)).*(taue.*taue2+etap.*taup.*(taue+2*taue2)).*(etay.*tauy.*taue.*taue2+etap.*taup.*(taue.*taue2+etay.*tauy.*(taue+taue2)))));
   
    for tt=1:rep
        % Initialization of weight matrix
        W_0=rand(N)*w_max*M/N*0.001; % random initial weights
        W_0=W_0-diag(diag(W_0));

        %Initialization of connectivity matrix
        W=W_0;

        % Initialization of inhibition
        inh=mean(sum(W,2))/(N-1);
        W_inh=ones(N)*inh;
        W_inh=W_inh-diag(diag(W_inh));

        % Initialization of heterosynaptic competition
        competition=zeros(N);
        max_competition=0.05*w_max; %upper boundary heterosynaptic competition
        
        % Initialize internal correlations terms:
        first_order = zeros(N);
        second_order = zeros(N);
        third_order = zeros(N);
        
        %% Starting the stimulations
        for ii=1:iterations

            W_eff=W-W_inh; % Effective excitatory connectivity
            W2=W_eff*W_eff;
            W3=W_eff*W2;

            % Calculate the average firing rate
            firing_rate=(eye(N)-eps_0*W_eff)^-1*b_0;
            D=diag(firing_rate);
            firing_rate_inv=1./firing_rate;
            D_inv=diag(firing_rate_inv);


            if sum((firing_rate)<0)>N/10
                error('firing rate error')
            end

            % Calculate the average change in connectivity weight for each order, at each iteration.
            % First order contribution:
            first_order = W_eff*D*M10+D*W_eff'*M01;

            % Second order contribution:
            second_order=W2*D*M20+D*W2'*M02+W_eff*D*W_eff'*M11+D_inv*repmat(diag(W_eff*D*W_eff'),1,N)*D*Ma1g1+D_inv*repmat(diag(W2),1,N)*D*Ma2g0+D_inv*W_eff.^2*D*Ma1b0g1;

            % Third order contribution:
            third_order=W3*D*M30+D*W3'*M03+W2*D*W_eff'*M21+W_eff*D*W2'*M12+D_inv*repmat(diag(W_eff*D*W2),1,N)*D*Ma2g1plusa1g2+D_inv*repmat(diag(W3),1,N)*D*Ma3g0+D_inv*W_eff.^2*D*W_eff'*Ma1b1g1;

            % Total change
            dW=real(first_order+second_order+third_order+extcorr.*M11);

            % Min bound of heterosynaptic competition
            if max_competition~=0
                mu=min([1./(max(max(abs(dW)))*nu)*w_max*0.0005 1./(max_competition)]);
            else
                mu=1./(max(max(abs(dW)))*nu)*w_max*0.0005;
            end

            % Adaptive numerical integration
            W=W+(mu*nu)*(dW);
            W=W-diag(diag(W));
            W=W.*(W>0).*(W<w_max)+(W>=w_max)*w_max; % hard bound

            % Heterosynaptic competition:
            sum_rows=sum(W,2)*ones(1,N);
            sum_column=ones(N,1)*sum(W);
            competition=-Psi*(sum_rows-W_max).*(W_max<sum_rows)-Psi*(sum_column-W_max).*(W_max<sum_column);

            W=W+mu*competition;
            W=W.*(W>0);
            W=W-diag(diag(W));
            
            max_competition=max(max(abs(competition))); % check the learning rate of the competition in order to ajust the numerical integration

            % Update inhibition:
            inh=mean(sum(W,2))/(N-1);
            W_inh=ones(N)*inh;
            W_inh=W_inh-diag(diag(W_inh));

        end
        Wt{kk,tt}=W; 
    end
end
save(savefile);
