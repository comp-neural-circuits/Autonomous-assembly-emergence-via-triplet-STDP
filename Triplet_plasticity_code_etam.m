%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is the basis of the manuscript:
% Lisandro Montangie, Christoph Miehl, Julijana Gjorgjieva (2020) PLoS
% Computational Biology
%
%
% The following code is used to calculate the evolution of the average
% synaptic efficacy in a recurrent network based on motif expansion up to
% third-order of the triplet STDP rule (Eqs. 44-48). A detailed
% derivation of the motif expansion is done in the manuscript.
% The code is based on code from Ravid Tannenbaum, Burak (2016) PLoS Comput Biol 12(8):e1005056
%
%
% With this code, Fig.6-8,11 and Fig. S2 from the manuscript can be
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
ETAM=[1:25]; % Vector of depression modulation parameters
inh_mult=1; % Here the detailed balance of excitation and inhibition can be changed
delta=0; % Here the balance of potentation and depression can be distrupted (ie zero-order motif non-zero)
iterations=200000; % Number of time-steps

N=48; % Number of neurons in the network
rep=100; % Number of repetitions
mu0=0.15; % Input firing rate (ms^(-1))
b_0=ones(N,1)*mu0;

% Heterosynaptic competition parameters
M=5; %Number of strong synapses
w_max=0.85/M; % keeps the maximum sum of each row and column synaptic efficacy constant so that spectral radius <1 and approximation is good
W_max=M*w_max;% maximum row and coloumn synaptic strength
Psi=0.1/(M*(w_max.^2)); % heterosynaptic competition rate
nu=Psi*0.0005; % learning rate

Wt=cell(length(ETAM),rep);

savefile='data_triplet_STDP';

%% Start simulation
for k= 1:numel(ETAM)
    etam=ETAM(k);
    
    % Calculation of motif coefficients (for details see Methods and Supplementary information 2)
    % Figure 5 shows all motifs used in this simulation (except the zero-order motif if delta is non-zero)
    M0=delta;

    deltac=1+delta/(Am*taum);

    % motif alpha=1, beta=0
    M10=2*pi.*Am*taum.*deltac.*(taup./((taue+taup).*(tauef+taup))+(tauy).^2./((taue+tauy).*(taup+tauy).*(tauef+tauy)));

    % motif alpha=0, beta=1
    M01=-2*pi.*Am*taum.*((etam.*taum)./((taue+etam.*taum).*(tauef+etam.*taum))-deltac.*(taup).^2./((taue+taup).*(taup+tauy).*(tauef+taup)));

    % motif alpha=1, beta=1
    M11=pi.*Am*taum.*(deltac.*(taup+tauef2)./((taue+taup).*(tauef+taup))+deltac.*2.*(taup).^4./((taup+tauy).*((taup).^2-taue.^2).*((taup).^2-tauef.^2)) ...
        +deltac.*(taue+taue2).^2./((taue+2.*taue).*(taue-taup).*(taue+tauy))+deltac.*taue2.^3./((taue+taue2).*(taup-tauef).*(tauy+tauef)) ...
        -(etam.*taum+tauef2)./((taue+etam.*taum).*(tauef+etam.*taum)));

    % motif alpha=2, beta=0
    M20=2*pi.*Am*taum.*deltac.*((taup).^3./((taue+taup).^2.*(tauef+taup).^2)+(tauy).^4./((taue+tauy).^2.*(taup+tauy).*(tauef+tauy).^2));

    % motif alpha=0, beta=2
    M02=-2*pi.*Am*taum.*((etam.*taum).^3./((taue+etam.*taum).^2.*(tauef+etam.*taum).^2)-deltac.*(taup).^4./((taue+taup).^2.*(taup+tauy).*(tauef+taup).^2));

    % motif alpha=2, beta=1
    M21=pi/2.*Am.*taum.*(deltac.*(13*(taup.*taue).^3+taue.^4.*(taue2+taup).^2+taue.^3.*(taue2+taup).*(2*taue2+taup).*(taue2+3.*taup)...
        +taup.*taue.^2.*taue2.*(3.*taue2+4.*taup).^2+2.*(taup.*taue2).^2.*taue.*(8*taue2+13*taup))./((taue+taue2).*(taue+2.*taue2).^2.*(taue+taup).^2.*(tauef+taup).^2)...
        + deltac.*(taue + taue2).*(tauy.*taue.^2.*taue2.*(taue.^3.*taue2.^2.*(taue+2.*taue2)+tauy.*taue.^2.*taue2.*(2*taue.^2+9*taue.*taue2+9*taue2.^2)+(tauy).^3.*(3*taue.^3+16*taue.^2.*taue2+26*taue.*taue2.^2 ...
        +13*taue2.^3)+(tauy).^2.*taue.*(taue.^3+10*taue.^2.*taue2+24.*taue.*taue2.^2+16.*taue2.^3))+taup.*taue.*(taue+2.*taue2).*(taue.^4.*taue2.^3+3.*tauy.*taue.^3.*taue2.^2.*(taue+2.*taue2)...
        +3.*(tauy).^2.*taue.^2.*taue2.*(taue.^2+5*taue.*taue2+5*taue2.^2)+(tauy).^4.*(3.*taue.^3+17*taue.^2.*taue2+28.*taue.*taue2.^2+14*taue2.^3)+(tauy).^3.*taue.*(taue.^3+12.*taue.^2.*taue2...
        +30*taue.*taue2.^2+20*taue2.^3))+(taup).^2.*(4*(tauy).^4.*(taue.^2+3.*taue.*taue2+2*taue2.^2).^2+taue.^4.*taue2.^2.*(taue.^2+3*taue.*taue2+3*taue2.^2)+tauy.*taue.^3.*taue2.*(2*taue.^3+11*taue.^2.*taue2...
        +21*taue.*taue2.^2 +14*taue2.^3)+(tauy).^2.*taue.^2.*(taue.^4+11*taue.^3.*taue2+38*taue.^2.*taue2.^2+54.*taue.*taue2.^3+27.*taue2.^4)+(tauy).^3.*taue.*(3*taue.^4+23*taue.^3.*taue2...
        +62.*taue.^2.*taue2.^2+70.*taue.*taue2.^3+28.*taue2.^4)))./((taup + tauy).*(taup+taue).*(tauy+taue).^2.*(taue+2.*taue2).^2.*(taue.*taue2+taup.*(taue+taue2)).*(taue.*taue2+tauy.*(taue+taue2)).^2)...
        -((taue+taue2).^3./(taue.*(taue+2.*taue2).^3.*(taue+etam.*taum))-taue2.^3./(taue.*(taue+2.*taue2).^2.*(tauef+etam.*taum))));
    
    % motif alpha=1, beta=2
    M12=-pi/2.*Am.*taum.*((13*(etam.*taum.*taue).^3+taue.^4.*(taue2+etam.*taum).^2+taue.^3.*(taue2+etam.*taum).*(2*taue2+etam.*taum).*(taue2+3.*etam.*taum)...
        +etam.*taum.*taue.^2.*taue2.*(3.*taue2+4.*etam.*taum).^2+2.*(etam.*taum.*taue2).^2.*taue.*(8*taue2+13*etam.*taum))./((taue+taue2).*(taue+2.*taue2).^2.*(taue+etam.*taum).^2.*(tauef+etam.*taum).^2)...
        -deltac.*(taue+taue2).*(tauy.*taue.^4.*taue2.^2.*(taue.*taue2.*(taue+2.*taue2)+tauy.*(taue.^2+3*taue.*taue2+3*taue2.^2))+taup.*taue.^3.*taue2.*(taue+2*taue2).*(taue.^2.*taue2.^2 ...
        +3*tauy.*taue.*taue2.*(taue+2*taue2)+(tauy).^2.*(2*taue.^2+7*taue.*taue2+7*taue2.^2))+(taup).^3.*taue.*(taue+2*taue2).*(taue.^2.*taue2.*(taue.^2+8.*taue.*taue2+8*taue2.^2) ...
        +(tauy).^2.*(3*taue.^3+17*taue.^2.*taue2+28*taue.*taue2.^2+14*taue2.^3)+tauy.*taue.*(taue.^3+12.*taue.^2.*taue2+30*taue.*taue2.^2+20*taue2.^3)) ...
        +(taup).^4.*(taue+taue2).*(4*(tauy).^2.*(taue+taue2).*(taue+2*taue2).^2+taue.^2.*taue2.*(3*taue.^2+13*taue.*taue2+13*taue2.^2)+tauy.*taue.*(3*taue.^3+20.*taue.^2.*taue2 ...
        +42*taue.*taue2.^2+28*taue2.^3))+(taup).^2.*taue.^2.*(taue.^2.*taue2.^2.*(2*taue.^2+9*taue.*taue2+9*taue2.^2)+3*tauy.*taue.*taue2.*(taue.^3+7*taue.^2.*taue2+15*taue.*taue2.^2+10*taue2.^3) ...
        +(tauy).^2.*(taue.^4+11*taue.^3.*taue2+38*taue.^2.*taue2.^2+54*taue.*taue2.^3+27*taue2.^4)))./((taup+tauy).*(taup+taue).^2.*(tauy+taue).*(taue+2*taue2).^2.*(taue.*taue2 ...
        +taup.*(taue+taue2)).^2.*(taue.*taue2+tauy.*(taue+taue2)))...
        -deltac.*((taue+taue2).^3./(taue.*(taue+2.*taue2).^3.*(taue+taup))-taue2.^3./(taue.*(taue+2.*taue2).^2.*(tauef+taup))));

    % motif alpha=3, beta=0
    M30=2*pi.*Am*taum.*deltac.*((taup).^5./((taue+taup).^3.*(tauef+taup).^3)+(tauy).^6./((taue+tauy).^3.*(taup+tauy).*(tauef+tauy).^3));

    % motif alpha=0, beta=3
    M03=-2*pi.*Am*taum.*((etam.*taum).^5./((taue+etam.*taum).^3.*(tauef+etam.*taum).^3)-deltac.*(taup).^6./((taue+taup).^3.*(taup+tauy).*(tauef+taup).^3));
    
    % motif alpha=2, gamma=0 from auto-covariance C_ii
    Ma2g0=2*pi.*Am*taum.*deltac.*((tauy).^3./((taue+tauy).^2.*(tauef+tauy).^2)); %the rate is considered in the sum for this terms

    % motif alpha=3, gamma=0 from auto-covariance C_ii
    Ma3g0=2*pi.*Am*taum.*deltac.*((tauy).^5./((taue+tauy).^3.*(tauef+tauy).^3)); %the rate is considered in the sum for this terms
    %the rate is considered in the sum for this terms

    % motif alpha=1, gamma=1 from auto-covariance C_ii
    Ma1g1=pi.*Am*taum.*deltac.*((tauy+tauef2)./((taue+tauy).*(tauef+tauy)));

    % motif alpha=2, gamma=1 and alpha=1, gamma=2 from auto-covariance C_ii
    Ma2g1plusa1g2=pi/2.*Am.*taum.*deltac.*((taue+taue2).*(2*taue.^4.*(taue2+tauy).^2+taue2.^3.*tauy.^2.*(-5*taup+21*tauy)...
        +2.*taue.*taue2.^2.*tauy.*(7*tauy.*(2*taue2+3*tauy)-taup.*(2*taue2+5.*tauy))...
        +taue.^3.*(taue2+tauy).*(4*taue2.^2+tauy.*(-taup+5*tauy)+taue2.*(-taup+13*tauy))...
        +taue.^2.*taue2.*(-taup.*(taue2.^2+6.*taue2.*tauy+6.*tauy.^2)+tauy.*(17*taue2.^2+42.*taue2.*tauy...
        +26.*tauy.^2))))./((taue+2*taue2).^2.*(taue+tauy).^2.*(taue2.*tauy+taue.*(taue2+tauy)).^2);

    % motif alpha=1, beta=0, gamma=1 from third order cumulant K_ij
    Ma1b0g1=4*pi^2.*Am.*taum.* deltac.*((taup).^3.*(taue+taue2).^2.*(3*tauy.*taue.*taue2+2*taup.*(taue*taue2+tauy.*(taue+2*taue2))))...
        ./((2*taup+taue).*(tauy.*taue+taup.*(tauy+taue)).*(taup.*tauy.*taue+tauy.*taue.*taue2...
        +taup.*(tauy+taue).*taue2).*(taue.*taue2+2.*taup.*(taue+taue2)).*(taue.*taue2+taup.*(taue+2*taue2)));

    % motif alpha=1, beta=1, gamma=1 from third order cumulant K_ij
    Ma1b1g1=(2/3)*pi^2*Am.*taum.*deltac.*(taue+taue2).^3.*(1/(taue.^2.*(2*taup+taue).*(-tauy+taue).*(taue+3.*taue2))+taue2.^3./(taue.^2.*(taue+taue2).*(2*taue+3*taue2).*(taue.*taue2 ...
        +2*taup.*(taue + taue2)).*(taue.*taue2-tauy.*(taue+taue2)))-(3.*taue2.^2.*(-2*taue.*taue2+tauy.*(taue+2*taue2)))./(taue.^2.*(-tauy+taue).*(taue ...
        +3*taue2).*(2*taue+3*taue2).*(-taue.*taue2+tauy.*(taue+taue2)).*(taue.*taue2+taup.*(taue+2*taue2)))-(2*tauy.*taue+7*tauy.*taue2+2*taue.*taue2) ...
        ./((taup-taue).*taue.*(2*tauy+taue).*(taue+3*taue2).*(2*taue+3*taue2).*(taue.*taue2+tauy.*(taue+2*taue2)))+(3.*(tauy).^4)./((tauy-taue)...
        .*taue.*(2*tauy+taue).*(tauy.*taue+taup.*(tauy+taue)).*(-taue.*taue2+tauy.*(taue+taue2)).*(taue.*taue2+tauy.*(taue+2*taue2))) ...
        -(taue2.^3.*(5*tauy.*taue+7.*tauy.*taue2+2*taue.*taue2))./(taue.*(taue+taue2).*(taue+3*taue2).*(2*taue+3*taue2).*(taue.*taue2-taup.*(taue+taue2)) ...
        .*(taue.*taue2+2*tauy.*(taue+taue2)).*(taue.*taue2+tauy.*(taue+2*taue2)))-(3*(tauy).^4.*taue2.^2)./((tauy-taue).*taue.*(-taue.*taue2+tauy.*(taue+taue2)).*(taue.*taue2 ...
        +2*tauy*(taue+taue2)).*(taue.*taue2+tauy.*(taue+2*taue2)).*(tauy.*taue.*taue2+taup.*(taue.*taue2+tauy.*(taue+taue2))))+(3*(taup).^5.*(3*tauy.*taue.*taue2 ...
        +2*taup*(taue.*taue2+tauy.*(taue+2*taue2))))./((taup-taue).*(2*taup+taue).*(tauy.*taue+taup.*(tauy+taue)).*(-taue.*taue2+taup.*(taue+taue2))...
        .*(taue.*taue2+2*taup.*(taue+taue2)).*(taue.*taue2+taup.*(taue+2*taue2)).*(tauy.*taue.*taue2+taup.*(taue.*taue2+tauy.*(taue+taue2)))));

    for tt=1:rep
        % Initialization of the connectivity matrix
        W_0=rand(N)*w_max*M/N*0.001; % random initial weights
        W_0=W_0-diag(diag(W_0));
        W=W_0;

        % Initialization of inhibition
        inh=mean(sum(W,2))/(N-1)*inh_mult;       
        W_inh=ones(N)*inh;
        W_inh=W_inh-diag(diag(W_inh));

        % Initialization of heterosynaptic competition
        competition=zeros(N);
        max_competition=0.05*w_max; %upper boundary on heterosynaptic competition

        % Initialize first, second and third order motif contributions
        first_order = zeros(N); 
        second_order = zeros(N); 
        third_order = zeros(N);

        %% Starting the simulations
        for ii=1:iterations

            W_eff=W-W_inh; % Effective excitatory connectivity
            W2=W_eff*W_eff;
            W3=W_eff*W2;

            % Calculate the average firing rate
            firing_rate=(eye(N)-eps_0*W_eff)^-1*b_0;                
            D=diag(firing_rate);
            firing_rate_inv=1./firing_rate;
            D_inv=diag(firing_rate_inv);

            % Calculate the average change in connectivity weight for each motif order, at each iteration.
            % First order contribution
            first_order = W_eff*D*M10+D*W_eff'*M01;

            % Second order contribution
            second_order = W2*D*M20+D*W2'*M02+W_eff*D*W_eff'*M11+D_inv*repmat(diag(W_eff*D*W_eff'),1,N)*D*Ma1g1+D_inv*repmat(diag(W2),1,N)*D*Ma2g0+D_inv*W_eff.^2*D*Ma1b0g1;

            % Third order contribution
            third_order = W3*D*M30+D*W3'*M03+W2*D*W_eff'*M21+W_eff*D*W2'*M12+D_inv*repmat(diag(W_eff*D*W2),1,N)*D*Ma2g1plusa1g2...
                +D_inv*repmat(diag(W3),1,N)*D*Ma3g0+D_inv*W_eff.^2*D*W_eff'*Ma1b1g1;

            % Total change
            dW=real(M0*(firing_rate*firing_rate')+first_order+second_order+third_order);

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
            inh=mean(sum(W,2))/(N-1)*inh_mult;
            W_inh=ones(N)*inh;  
            W_inh=W_inh-diag(diag(W_inh));

        end
        Wt{k,tt}=W; 
    end       
end
save(savefile);