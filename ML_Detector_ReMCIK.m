% ============Maximum Likelihood Detector
function [BB,MM] = ML_Detector_ReMCIK(avSNR, M,p1,PwrSC,index_all,y,h,N,jj,ref_symmm,ref_sym)

dis_all = zeros(2.^p1,M); % all possible realizations, index & symbol
for mm = 1:M    
    % Average Energy Per QAM Symbol
    sym_norm = ref_symmm(mm,:);  
    sym_m = sqrt(PwrSC)*sym_norm;
    dis_m = zeros(2^p1,1);
    for bb = 1:2^p1
        pp = index_all(bb,:)+1;
        sym_b = zeros(N,1);
        sym_b(pp) = sym_m;
        % distance for all possible symbols      
            tmp = (y(:,jj)-avSNR*h(:,jj).*sym_b)'*(y(:,jj)-avSNR*h(:,jj).*sym_b); 
            %tmp = (y(:,jj)-h(:,jj).*sym_b)'*(y(:,jj)-h(:,jj).*sym_b); 
        dis_m(bb) = tmp;
    end
    dis_all(:,mm) = dis_m;
end
% find minimun row (symbol)
[tmp0, I] = min(dis_all);
% find minimum column (index)
[~, min_col] = min(tmp0);
% index demodulate
BB = I(min_col);
% M-ary demodulate
MM = ref_sym(min_col);