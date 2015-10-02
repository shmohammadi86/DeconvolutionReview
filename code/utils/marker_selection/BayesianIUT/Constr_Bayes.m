function r=Constr_Bayes(DATA,nrreplics,target,updown,nrruns)
%CONSTR_BAYES    Bayesian procedure for inequality constraints
%   r = CONSTR_BAYES(DATA,nrreplics,target,updown,nrruns) gives a Bayesian 
%   evaluation of the constrained hypothesis 
%   H1: mu_i>mu_j for j=1...J (or, mu_i<mu_j for j=1...J) versus
%   H2: not H1
%   INPUT: DATA: gene expression values in a probe by array matrix; arrays are
%   supposed to be grouped (replications are found in succeeding columns)
%       NRREPLICS: number of replications per condition
%       TARGET: condition supposed to be overexpressed / underexpressed
%       UPDOWN: binary indicator for up- (1) or down-regulation (-1)
%       of the target
%       NRRUNS: number of Gibbs samples
%   OUTPUT: Vector of Bayes factors, one value per gene
%   TOOLBOX: Uses the statistical toolbox
%Using this file implies that you agree with the license (see license.txt)

%Author: K. Van Deun, Department of Psychology, Catholic University of
%Leuven (Belgium)


[nrprobesets nrarrays]=size(DATA);
nrconditions=size(nrreplics,2);

%Preliminary: Calculation of condition-specific averages for all probesets
lower=1;
AV1=[];
AV2=[];
for i=1:nrconditions
    upper=lower+nrreplics(i)-1;
    v=ones(nrreplics(i),1);
    average=DATA(:,lower:upper)*v/(v'*v);
    AV1(:,lower:upper)=average*v';
    AV2=[AV2 average];
    lower=upper+1;
end;

%1. Initial values
sigma2_0=sum((AV1-DATA).^2,2)/nrarrays;
SIGMA2_SAMPLED=sigma2_0*ones(1,nrconditions);

%2. Gibbs Sampler
indicsum=zeros(nrprobesets,1);
NRREPLICS=ones(nrprobesets,1)*nrreplics;
for i=1:nrruns
    SEM=sqrt(SIGMA2_SAMPLED./NRREPLICS);
    MU_SAMPLED=AV2+(SEM.*randn(nrprobesets,nrconditions));
    lower=1;
    for j=1:nrconditions
        upper=lower+nrreplics(j)-1;
        MU_SAMPLED2(:,lower:upper)=MU_SAMPLED(:,j)*ones(1,nrreplics(j));
        lower=upper+1;
    end;
    SIGNDIFF=updown*sign(MU_SAMPLED(:,target)*ones(1,nrconditions)-MU_SAMPLED);
    indic=sum(SIGNDIFF,2)==(nrconditions-1);
    indicsum=indicsum+indic;
    df=nrarrays-2;
    nscale=sum((MU_SAMPLED2-DATA).^2,2);
    chisqvals=2.*randg(df./2,nrprobesets,1);
    sigma2_sampled=nscale.*chisqvals.^(-1);
    SIGMA2_SAMPLED=sigma2_sampled*ones(1,nrconditions);
end;
f1=indicsum/nrruns;
f2=1-f1;
BF=nrconditions*f1./((nrconditions/(nrconditions-1))*f2);
r=BF;