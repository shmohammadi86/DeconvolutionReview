function [ paths ] = setDatasetPaths( datasetName )
    switch(datasetName)                  
        case 'LiverBrainLung'
            paths.expression.mixture='./input/LiverBrainLung_GSE19830/mix.txt';
            paths.expression.signature='./input/LiverBrainLung_GSE19830/sig.txt';
            paths.expression.pure='./input/LiverBrainLung_GSE19830/pure.txt';

            paths.annotations.references='./input/LiverBrainLung_GSE19830/reference_annotations.txt';
            paths.annotations.features='./input/LiverBrainLung_GSE19830/feature_annotations.txt';
            paths.annotations.pure='./input/LiverBrainLung_GSE19830/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/LiverBrainLung_GSE19830/coef.txt';
            
        case 'BreastBlood'
            paths.expression.mixture='./input/BreastBlood_GSE29832/mix.txt';
            paths.expression.signature='./input/BreastBlood_GSE29832/sig.txt';
            paths.expression.pure='./input/BreastBlood_GSE29832/pure.txt';
            
            paths.annotations.references='./input/BreastBlood_GSE29832/reference_annotations.txt';
            paths.annotations.features='./input/BreastBlood_GSE29832/feature_annotations.txt';
            paths.annotations.pure='./input/BreastBlood_GSE29832/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/BreastBlood_GSE29832/coef.txt';

        case 'CellLines'
            paths.expression.mixture='./input/CellLines_GSE11058/mix.txt';
            paths.expression.signature='./input/CellLines_GSE11058/sig.txt';
            paths.expression.pure='./input/CellLines_GSE11058/pure.txt';

            paths.annotations.references='./input/CellLines_GSE11058/reference_annotations.txt';
            paths.annotations.features='./input/CellLines_GSE11058/feature_annotations.txt';
            paths.annotations.pure='./input/CellLines_GSE11058/pure_annotations.txt';            
            paths.annotations.samples='';

            paths.coeff='./input/CellLines_GSE11058/coef.txt';

        case 'CellLines_minimal'
            paths.expression.mixture='./input/CellLines_GSE11058/mix_minimal.txt';
            paths.expression.signature='./input/CellLines_GSE11058/sig_minimal.txt';
            paths.expression.pure='./input/CellLines_GSE11058/pure_minimal.txt';

            paths.annotations.references='./input/CellLines_GSE11058/reference_annotations.txt';
            paths.annotations.features='./input/CellLines_GSE11058/feature_annotations.txt';
            paths.annotations.pure='./input/CellLines_GSE11058/pure_annotations.txt';            
            paths.annotations.samples='';

            paths.coeff='./input/CellLines_GSE11058/coef.txt';            

        case 'RatBrain'
            paths.expression.mixture='./input/RatBrain_GSE19380/mix.txt';
            paths.expression.signature='./input/RatBrain_GSE19380/sig.txt';
            paths.expression.pure='./input/RatBrain_GSE19380/pure.txt';

            paths.annotations.references='./input/RatBrain_GSE19380/reference_annotations.txt';
            paths.annotations.features='./input/RatBrain_GSE19380/feature_annotations.txt';
            paths.annotations.pure='./input/RatBrain_GSE19380/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/RatBrain_GSE19380/coef.txt';
            
        case 'MAQC'
            paths.expression.mixture='./input/MAQC_GSE5350/mix_minimal.txt';
            paths.expression.signature='./input/MAQC_GSE5350/sig.txt';
            paths.expression.pure='./input/MAQC_GSE5350/pure.txt';

            paths.annotations.references='./input/MAQC_GSE5350/reference_annotations.txt';
            paths.annotations.features='./input/MAQC_GSE5350/feature_annotations.txt';
            paths.annotations.pure='./input/MAQC_GSE5350/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/MAQC_GSE5350/coef_minimal.txt';

        case 'Retina'
            paths.expression.mixture='./input/Retina_GSE33076/mix.txt';
            paths.expression.signature='./input/Retina_GSE33076/sig.txt';
            paths.expression.pure='./input/Retina_GSE33076/pure.txt';

            paths.annotations.references='./input/Retina_GSE33076/reference_annotations.txt';
            paths.annotations.features='./input/Retina_GSE33076/feature_annotations.txt';
            paths.annotations.pure='./input/Retina_GSE33076/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/Retina_GSE33076/coef.txt';
            
        case 'PERT_Cultured'
            paths.expression.mixture='./input/PERT_Cultured_GSE16589/mix.txt';
            paths.expression.signature='./input/PERT_Cultured_GSE16589/sig.txt';
            paths.expression.pure='./input/PERT_Cultured_GSE16589/pure.txt';
            
            paths.annotations.references='./input/PERT_Cultured_GSE16589/reference_annotations.txt';
            paths.annotations.features='./input/PERT_Cultured_GSE16589/feature_annotations.txt';
            paths.annotations.pure='./input/PERT_Cultured_GSE16589/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/PERT_Cultured_GSE16589/coef.txt';
            
        case 'PERT_Uncultured'
            paths.expression.mixture='./input/PERT_uncultured_GSE40830/mix.txt';
            paths.expression.signature='./input/PERT_uncultured_GSE40830/sig.txt';
            paths.expression.pure='./input/PERT_uncultured_GSE40830/pure.txt';
            
            paths.annotations.references='./input/PERT_uncultured_GSE40830/reference_annotations.txt';
            paths.annotations.features='./input/PERT_uncultured_GSE40830/feature_annotations.txt';
            paths.annotations.pure='./input/PERT_uncultured_GSE40830/pure_annotations.txt';
            paths.annotations.samples='';

            paths.coeff='./input/PERT_uncultured_GSE40830/coef.txt';


%         case 'Abbas+tumor+noise'
%             paths.expression.mixture='./input/Abbas+tumor+noise/mix.txt';
%             paths.expression.signature='./input/Abbas+tumor+noise/sig.txt';
%             paths.coeff= './input/Abbas+tumor+noise/coef.txt';
%             paths.annotations.references='./input/Abbas+tumor+noise/reference_annotations.txt';
%             paths.annotations.features='';
%             paths.annotations.samples='./input/Abbas+tumor+noise/sample_annotations.txt';
% 
% 
%         case 'PBMC_LM22'
%             paths.expression.mixture='./input/PBMC_Alizadeh/mix.txt';
%             paths.expression.signature='./input/PBMC_Alizadeh/sig_LM22.txt';
%             paths.coeff='./input/PBMC_Alizadeh/coef.txt';
%             paths.annotations.references='./input/PBMC_Alizadeh/reference_annotations.txt';
%             paths.annotations.features='';
%             paths.annotations.samples='';
            
        otherwise
            error('Unknown dataset %s\n', datasetName);
            paths.expression.mixture='';
            paths.expression.signature='';
            paths.coeff='';
            paths.annotations.references='';
            paths.annotations.features='';
            paths.annotations.samples='';
    end
end

