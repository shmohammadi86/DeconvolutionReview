function [yesno pval]=IUT(DATA,nrreplics,target,updown,alpha,mc)
%IUT    Calculates the Intersection Union Test in an ANOVA setting
%   [yesno pval] = IUT(DATA,nrreplics,target,updown,alpha,mc)
%   INPUT: DATA: gene expression values in a probe by array matrix; arrays are
%   supposed to be grouped (replications are found in succeeding columns)
%       NRREPLICS: number of replications per condition
%       TARGET: condition supposed to be overexpressed / underexpressed
%       UPDOWN: binary indicator for up- (1) or down-regulation (-1)
%       of the target
%       ALPHA: significance level
%       MC: Multiple comparisons corrections (1) or not (0)
%   OUTPUT: yesno: Vector of IUT binary outcomes (one value per gene)
%        pval:  p-value of the least significant TARGET-other condition pair
%   TOOLBOX: uses the statistical toolbox
%   References:
%   Berger, R.L. (1982). Multiparameter hypothesis testing and acceptance
%   sampling. Technometrics, 24, 295-300.
%Using this file implies that you agree with the license (see license.txt)

%Author: K. Van Deun, Department of Psychology, Catholic University of
%Leuven (Belgium)

[nrprobesets nrarrays]=size(DATA);
nrconditions=length(nrreplics);

%Calculation of condition-specific averages and variances for all probesets
lower=1;
AV=[];
SSWITHIN=0;
for i=1:nrconditions
    upper=lower+nrreplics(i)-1;
    v=ones(nrreplics(i),1);
    DATA_i=DATA(:,lower:upper);
    average=DATA_i*v/(v'*v);
    AV=[AV average];
    sswithin=((DATA_i-average*v').^2)*v;
    SSWITHIN=SSWITHIN+sswithin;
    lower=upper+1;
end;
%Calculation of t-statistics with s² from ANOVA context
MTARGET=AV(:,target)*ones(1,nrconditions-1);
MOTHERS=AV;
MOTHERS(:,target)=[];
invN_target=1/nrreplics(target);
invN_OTHERS=nrreplics.^(-1);
invN_OTHERS(target)=[];
SE2=(SSWITHIN/(nrarrays-nrconditions))*(invN_target+invN_OTHERS);
SE=SE2.^(.5);
TMATRIX=(MTARGET-MOTHERS)./SE;

%Critical t-value using multiple comparison correction for testing if mc=1
%multiple probes (Sidak ~= Bonferonni)
if mc==1
    alpha=1-(1-alpha)^(1/nrprobesets);%Sidak
end;
critT=tinv(1-alpha,nrarrays-nrreplics(target));

%SUMMARY
SUMSTAT=sum(updown*TMATRIX>critT,2);
minT=min(updown*TMATRIX,[],2);
pval = 1-tcdf(minT,nrarrays-nrreplics(target));
yesno=SUMSTAT==nrconditions-1;