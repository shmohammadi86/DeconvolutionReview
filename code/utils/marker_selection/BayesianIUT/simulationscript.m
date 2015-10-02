%SIMULATIONSCRIPT    Script for a Monte Carlo Simulation 
%      SIMULATIONSCRIPT
%           Evaluates and compares the procedures to test for 
%           tissue-selectivity (IUT & Bayesian method).
%Using this file implies that you agree with the license (see license.txt)

%Author: K. Van Deun, Department of Psychology, Catholic University of
%Leuven (Belgium)

tic

%The data are generated in maximal agreement with the microarray tissue
%data (GSE9954).
NRREPLICS=[3 3 4 3 3 3 3 3 4 3 5 3 3 3 3 3 3 3 3 3 3 3];
nrconditions=length(NRREPLICS);
stddev=sqrt(0.04);
av_target=7.3;
%The target tissue (first tissue) is supposed to be upregulated
TARGET=1;
UPDOWN=1;

%Factors of the simulation study:
%1. Nr of tissues with higher expression than the target tissue
% Two cases for one tissue: neutral (mu=mu_target) and not supporting
% (mu>mu_target)
nrnotsupport=[0 1 1 21];
%2. Effect size (Cohen's d: (diff means)/sqrt(pooled variance)).
delta=[0.5 2];

for i=1:length(nrnotsupport)
    for j=1:length(delta)
        
        
        %GENERATION OF THE DATA
        %The data are generated in maximal agreement with the microarray tissue
        %data (GSE9954)
        if i==3
            av_other1=av_target;%NEUTRAL case
            fprintf('\n NEUTRAL CASE \n')
        else
            fprintf('\n Number of tissues not supporting tissue selectivity: %3.0f \n',nrnotsupport(i))
            av_other1=av_target+delta(j)*stddev;
        end;
        av_other2=av_target-delta(j)*stddev;
        av=[av_target*ones(1,NRREPLICS(1)) av_other1*ones(1,sum(NRREPLICS(2:nrnotsupport(i)+1))) ...
            av_other2*ones(1,sum(NRREPLICS(nrnotsupport(i)+2:nrconditions)))];
        n_MC=5000;%Resampling size for Monte Carlo simulation
        AV=ones(n_MC,1)*av;
        DATA=AV+(stddev*randn(n_MC,sum(NRREPLICS)));
        
        %BAYESIAN ANALYSIS
        NRRUNS=5000;
        BF=Constr_Bayes(DATA,NRREPLICS,TARGET,UPDOWN,NRRUNS);
        k=find(BF>1);
        fprintf('Effect size: %3.3f \t Proportion BF > 1: %6.5f \n',delta(j),length(k)/n_MC)
        k=find(BF>32);
        fprintf('Effect size: %3.3f \t Proportion BF > 32: %6.5f \n',delta(j),length(k)/n_MC)
        
        %INTERSECTION UNION TEST (R.L. Berger, 1987, Technometrics)
        ALPHA=.05;
        MC=0;%no correction for multiple testing
        [BERGER PVAL]=IUT(DATA,NRREPLICS,TARGET,UPDOWN,ALPHA,MC);
        k=find(BERGER(:,1)>0);
        fprintf('Effect size: %3.3f \t Proportion Berger significant: %6.5f \n',delta(j),length(k)/n_MC)
    end;
end;
toc