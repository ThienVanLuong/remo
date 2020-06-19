%============ Repeated MCIK-OFDM (ReMO) ===================
% Created by Thien Van Luong, Queen's University Belfast, UK, in June 2017
% Email: tluong01@qub.ac.uk
% Publication: T. V. Luong, Y. Ko and J. Choi, “Repeated MCIK-OFDM with enhanced transmit diversity under CSI uncertainty“,
% IEEE Trans. Wireless Commun., Mar 2018.

clear
M=16;
N=4;
K=3;

CSI = 1;
var=0.05;
ro=1;
Plot_type = 2;

Mary=2; % 1 PSK, 2 QAM

if(M==8)
    QAM = (5*M-4)./6;
else
    QAM = (2/3)*(M-1);
end

if(M==2)
    xi=1;
else
    xi=2;
end

tic
%% ======================= MCIK Parameters ================================
iter = 4;  % Iterations
nSymPerFrame = 1e4;
% Number of symbol per frame(1 OFDM symbol)
EbN0dB = 0:5:35;
EbN0 = 10.^(EbN0dB/10);
% Es/N0 parameter
PwrSC = N/K; %N/K; % Average Tx power per active sub-carrier
%PwrSC=1;
bps = log2(M); % bits per symbol
EsN0dB = EbN0dB; % + 10*log10(bps);%+10*log10(1/PwrSC);
EsN0 = 10.^(EsN0dB/10);
c = 2^floor(log2(nchoosek(N,K)))-0; % Effective Carrier Combinations
p1 = floor(log2(nchoosek(N,K)))-0;  % index bit length per cluster
p2 = bps; % information bit length per cluster
sigma = sqrt(1./EsN0);
%%% M-PSK
sym_test=zeros(1,M);
for pp=1:M
    if(Mary==1)
        sym_test(pp)=pskmod(M-pp,M,ro*pi./M,'gray');
    else
        sym_test(pp)=qammod(M-pp,M,0*pi./M,'gray');
    end
end
ref_sym = sym_test.';
if(Mary==1)
    ref_symmm = ref_sym.*(1./abs(ref_sym));
else
    ref_symmm = ref_sym.*(1/sqrt(QAM));
end
%  index_all = Combin_Md(N,K);
if(K==2&&N==4)
    index_all = [1 0;2 0;3 1;3 2];
else
    index_all = Combin_Md(N,K);
end
index_allz=index_all+1;
MCIK_LUT=zeros(c,N);
for row=1:c
    for col=1:K
        MCIK_LUT(row,index_allz(row,col))=1;
    end
end
M1=2.^p1;
PP=zeros(M1,M1);
for i1=1:M1
    for j1=1:M1
        PP(i1,j1)=sum(abs(MCIK_LUT(i1,:)-MCIK_LUT(j1,:)));
    end
end
t1=sum(sum(PP==2))./M1;
t2=sum(sum(PP==4))./M1;
t3=t1+t2;

%% ==================== Loop for SNR =========================
PEP = zeros(1,size(sigma,2)); % index symbol error IEP
OFDM_SER = zeros(1,size(sigma,2)); % ofdm symbol error
Total_SER = zeros(1,size(sigma,2));
BER=zeros(1,size(sigma,2));
BER1=zeros(1,size(sigma,2)); % 
BER2=zeros(1,size(sigma,2));

for s1 = 1:size(sigma,2)
    fprintf('== EbN0(dB) is %g == \n',EbN0dB(s1))
    %% ==================== Loop for iteration =======================
    symerr_mcik = zeros(1,iter);
    symerr_ofdm = zeros(1,iter);
    symerr_iter= zeros(1,iter);
    BER_iter= zeros(1,iter);
    BER_iter_1= zeros(1,iter); % index bit error rate
    BER_iter_2= zeros(1,iter); % symbol bit error rate
    for s2 = 1:iter
        fprintf('== EbN0(dB) is %g and iteration is %g == \n',EbN0dB(s1),s2)
        %% ===================== Bit generator =========================
        % bit = (index bit + M-ary bps) * symbols in OFDM frame
        bit = randi([0 1],1,(p1+p2)*nSymPerFrame);
        % bit split - reshape bit stream (p1+p2)
        bit2 = reshape(bit.',p1+p2,nSymPerFrame).';
        %% ================= Index selector =========================
        % information bits (p2)
        info_bit = bit2(:,p1+1:end);
        % mapping bit to QAM symbol
        info_dec_i = bi2de(info_bit);
        if(Mary==1)
            sym = pskmod(info_dec_i,M,ro*pi./M,'gray');
            sym_norm = sym.*(1./abs(sym));
        else
            sym = qammod(info_dec_i,M,0*pi./M,'gray');
            sym_norm = sym.*(1/sqrt(QAM));
        end
        % index bits (p1)
        index_bit = bit2(:,1:p1);
        % index symbol ( bit to decimal ), select indices from combinatorial method
        index_sym = BitoDe(index_bit);
        % Power reallocation
        sym_tx = sym_norm.*sqrt(PwrSC);
        % transmitted OFDM symbols
        tx_sym = zeros(N,nSymPerFrame);
        for kk = 1:nSymPerFrame
            % kk-th index symbol for cluster
            kk_index = index_sym(kk)+1;
            % select combination
            indices = index_all(kk_index,:)+1;
            tx_sym(indices,kk) = sym_tx(kk,:);
        end
        %%%%%====================== Fixed expsilon^2 ===================
        if(CSI==1)
            eps=0;
        elseif(CSI==2)
            eps=var;
        else
            eps=1./(1+EsN0(s1));
        end
        
        noise = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)));
        h = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)))*sqrt(1-eps);
        e=sqrt(eps)./sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)));
        h1=h+e;
        y = sqrt(EsN0(s1))*h1.*tx_sym+noise;
        avSNR=sqrt(EsN0(s1));
        
        %% ================== ML detect ====================
        index_sym_de = zeros(1,nSymPerFrame);
        indices_de = zeros(nSymPerFrame,K);
        re_sym = zeros(nSymPerFrame,1);
        OutputData = zeros(N,nSymPerFrame);
        for jj=1:nSymPerFrame
            [BB,MM] = ML_Detector_ReMCIK(avSNR,M,p1,PwrSC,index_all,y,h,N,jj,ref_symmm,ref_sym);
            index_sym_de(jj) = BB-1;
            re_sym(jj,:) = MM;
        end
        %% =================error rate computation====================
        % ofdm symbol error
        ofdm_symerr = sum(sum(sym~=re_sym)); %
        % index symbol error
        ind_symerr = sum(index_sym~=index_sym_de);
        % index symbol to bit, index bit error
        index_bit_de = DetoBit(index_sym_de,p1);
        index_bit_err=sum(sum(index_bit~=index_bit_de));
        if(Mary==1)
            info_de_re=pskdemod(re_sym,M,ro*pi./M,'gray');
        else
            info_de_re=qamdemod(re_sym,M,0*pi./M,'gray');
        end
        info_bit_re= zeros(nSymPerFrame,bps);
        for kk=1:1
            info_bit_re(:,(kk-1)*bps+1:kk*bps)=de2bi(info_de_re(:,kk),bps);
        end
        info_bit_err=sum(sum(info_bit~=info_bit_re));
        %% ===========symbol & bit error rate  1 iteration==========
        % MCIK sym error
        symerr_mcik(s2) = ind_symerr./nSymPerFrame;
        % OFDM sym error
        symerr_ofdm(s2) = ofdm_symerr./nSymPerFrame;
        % symbol error rate
        symerr_iter(s2) = (ind_symerr+ofdm_symerr)./(2*nSymPerFrame);
        
        BER_iter(s2)=(info_bit_err+index_bit_err)./((p1+p2)*nSymPerFrame);
        BER_iter_1(s2) = index_bit_err./p1./nSymPerFrame;
        BER_iter_2(s2) = info_bit_err./p2./nSymPerFrame;
        
        
    end
    %% =============average symbol/bit error rate================
    PEP(s1) = sum(symerr_mcik)/iter; % IEP
    OFDM_SER(s1) = sum(symerr_ofdm)/iter;
    Total_SER(s1) = sum(symerr_iter)/iter;
    BER(s1)= sum(BER_iter)./iter;
    BER1(s1)= sum(BER_iter_1)./iter;
    BER2(s1)= sum(BER_iter_2)./iter;
end
fprintf('¬¬ N = %g / K = %g / M = %g / ¬¬ \n',N,K,M)

%% display plot figure
if Plot_type == 1
    plot = PEP;
else
    plot = Total_SER;
end


%% Theoretical Bound
SNR_av=EsN0*PwrSC;
PEP_theo=zeros(size(EbN0dB));
SEP_BMO=zeros(size(EbN0dB));
SEP_asym=zeros(size(EbN0dB));
PM=zeros(size(EbN0dB));
BER_theo=zeros(size(EbN0dB));
MM=(sin(pi./M)).^2;
if(K==2)
    xu=1;
else
    xu=0;
end
% t1=2.^p1-1;
%t1=2;
for i=1:length(EbN0dB)
    snr=SNR_av(i);
    Es = EsN0(i);
    if(CSI==1)
        eps=0;
    elseif(CSI==2)
        eps=var;
        a=1-eps;
    else
        eps=1./(1+EsN0(i));
    end
    PM(i)=xi*(1./(1+(1-eps)*MM*snr./(1+eps*snr./1)).^K./12+1./(1+4./3*MM*(1-eps)*snr./(1+eps*snr./1)).^K./4);
    PEP_theo(i)=(t1./12)*(1./(1+(1-eps)*snr./(4*(1+eps*snr./2))).^2+3./(1+(1-eps)*snr./(3*(1+eps*snr./2))).^2);
    SEP_BMO(i)=(PEP_theo(i)+PM(i))./2;
    BER_theo(i)=((p1./2+0)*PEP_theo(i)+PM(i))./(p1+p2);
    
    
    t2=xu./8./MM.^2;
    if(CSI==1)
        SEP_asym(i)=43./24.*(t1+t2)./snr.^2;
    elseif(CSI==2)
        SEP_asym(i)=(t1./24).*((1+a./2./eps).^(-2)+3.*(1+2*a./3./eps).^(-2))+((1+a*MM./eps).^(-K)+3.*(1+4*a*MM./3./eps).^(-K))./12;
    else
        SEP_asym(i)=43./24.*(t1.*(K./N+0.5).^2+t2.*(K./N+1).^2)./Es.^2;
    end
    
end


figure (23)
%semilogy(EbN0dB,plot,'b O-','LineWidth',1.5,'MarkerSize',10)
%hold on
% semilogy(EbN0dB,PEP_theo,'k :','LineWidth',1.5,'MarkerSize',10)

semilogy(EbN0dB,BER,'b +-','LineWidth',1.5,'MarkerSize',10)
hold on
%  semilogy(EbN0dB,BER_theo,'k :','LineWidth',1.5,'MarkerSize',10)
%  hold on
%   semilogy(EbN0dB,SEP_BMO,'k :','LineWidth',1.5,'MarkerSize',10)
%  hold on
%  semilogy(EbN0dB,SEP_asym,'k--','LineWidth',1.5)
%  hold on
% semilogy(EbN0dB,PEP,'b *-')
% hold on
% semilogy(EbN0dB,PEP_theo,'b--o')
axis([0 40 10^-5 10^0])
grid on
hold on
title('')
xlabel('Es/No (dB)')
if Plot_type == 1
    ylabel('Average IEP')
else
    ylabel('Average SEP')
end
toc