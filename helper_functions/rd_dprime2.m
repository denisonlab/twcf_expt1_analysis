function [dprime, criterion] = rd_dprime2(nh,nfa,nsignal,nnoise)

% [dprime, criterion] = rd_dprime2(nh,nfa,nsignal,nnoise)
%
% nh = number of hit trials
% nfa = number of false alarm trials
% nsignal = number of signal trials
% nnoise = number of noise trials

% adjust for ceiling or floor values
% nh(nh==nsignal) = nsignal(nh==nsignal)-1;
% nh(nh==0) = 1;
% nfa(nfa==nnoise) = nnoise(nfa==nnoise)-1;
% nfa(nfa==0) = 1;

% proportions
h = nh./nsignal;
fa = nfa./nnoise;

% % loglinear adjustment (Stanislaw & Todorov 1999) 
if h==0 || h==1 || fa==0 || fa==1
    nh = nh + 0.5; 
    nfa = nfa + 0.5; 
    nnoise = nnoise + 1; 
    nsignal = nsignal + 1; 
end

% adjustment = 0.01; 
% % adjust for ceiling or floor values 
% if h==1
%     h = 1-adjustment; 
% end
% if fa==1
%     fa = adjustment; 
% end 
% 
% if h==0
%     h = adjustment; 
% end
% if fa==0
%     h = 1-adjustment; 
% end

% dprime
zh = norminv(h,0,1); zfa = norminv(fa,0,1);

dprime = zh - zfa;
criterion = -0.5*(zh+zfa);
