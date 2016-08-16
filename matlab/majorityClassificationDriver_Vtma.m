clc
clear all
close all

%%
% Add my code folders.
addpath('./functionMfiles/dataManagement')

%%
[isOctave,matType,matExt] = checkIsOctave();
  
%%
% Calibrant Details.
% calibrantMZ = [1296.685,1570.677,2147.199,2932.588];
% calibrantName = {'Angiotensin I','[Glu]-Fibrinopeptide B', ... 
%     'Dynorphin A','ACTH fragment (1-24)'};


%%
% Proccessing parameters
curDataParams = struct();
curDataParams.binSize = 0.25;
curDataParams.wiggle = true;
curDataParams.wiggle_den = 3;
curDataParams.ppmNormTol = true;
curDataParams.normTol = 10;

vFolNams = {'Vtma'};        

vDatTypes = {'Binary','Intensity','LogIntensity','SN','Area'};
vVRmethods = {'None' 'PCA' 'CCA'};
% vVRmethods = {'None' 'PCA' 'CCA' 'altCCA'};

vUseAnnot = [false true];

vNcal = 0;
% vNcal = [0 3 4];
v_minNspec_r = 1;
% v_minNspec_r = [1 100];
v_minNspec_bin = 1;
% v_minNspec_bin = [1 100];

fid = fopen('.\output\Vtma_majority_classification.txt','wt');

% Write Header
fprintf(fid,'Me,Mloo,N,Dataset,VRmethod,DAmethod,DataType,Normalisation,includeEmptyValues,Smooth,CancerAnnotation,minNcal,minNspecPerCore,minNspecPerBin,nComponents\n');

%% Dataset
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};

%% Data Type
for datType_idx = 1:length(vDatTypes) 
curDataParams.dataType = vDatTypes{datType_idx};
disp([curDataParams.dataType ' Data'])

%% Annotated Cancer Spectra
for use_cancer_annot = vUseAnnot
curDataParams.use_cancer_annot = use_cancer_annot;
if curDataParams.use_cancer_annot
    disp(['  ' 'Annotated Spectra Only'])
else
    disp(['  ' 'All Spectra'])
end    

%% Calibrant QA filter
for nCal_tol = vNcal
curDataParams.nCal_tol = nCal_tol;

%% Number of spectra per region filter
for nSpec_r = v_minNspec_r
curDataParams.nSpec_r_tol = nSpec_r;

%% Number of spectra per bin filter
for nSpec_bin = v_minNspec_bin
curDataParams.nSpec_bin_tol = nSpec_bin;

switch curDataParams.dataType
    case 'Binary'
        vSmoothParams = [0 0.15 0.25];
        vNorm = [false];
        vUseEmptyVals = [true];
        
    case {'LogIntensity','Intensity','SN','Area'}
        vSmoothParams = [0];
        vNorm = [false true];
        vUseEmptyVals = [true false];
        
end


%% Binary Smoothing
for s = vSmoothParams
curDataParams.smoothParam = s;

%% Normalisation
for n = vNorm
curDataParams.normalisation = n;

%% Include empty values when taking averages?
for includeEmptyVals = vUseEmptyVals
curDataParams.includeEmptyVals = includeEmptyVals;

%% Variable Reduction Method
for vr_method_idx = 1:length(vVRmethods)
vr_method_cur = vVRmethods{vr_method_idx};
disp(['  ' '  ' vr_method_cur])
     
switch vr_method_cur
    case 'None'
        vDAmethods = {'NB','DWD'};
    case {'PCA','CCA','altCCA'}
        vDAmethods = {'NB','LDA','DWD'};
end

%% Classification Method
for da_method_idx = 1:length(vDAmethods)
da_method_cur = vDAmethods{da_method_idx};
disp(['  ' '  ' '  ' da_method_cur])
    
switch vr_method_cur
    case 'None'
        curDataParams.classificationMethod = da_method_cur;
        vNcomponents = -1;
        
    case 'PCA'
        curDataParams.classificationMethod = ['pca' da_method_cur];
        curDataParams.restrict_p_suit = true;
        [mFileNam,~] = matFileNamSelect('pPCA',curDataParams);
        load([mFileNam matExt],'veigval')
        vNcomponents = 1:length(veigval);
        clear veigval
        
    case 'CCA'
        curDataParams.classificationMethod = ['cca2' da_method_cur];
        vNcomponents = 1:30;
        
    case 'altCCA'
        curDataParams.classificationMethod = ['cca1' da_method_cur];
        vNcomponents = 1:30;

end

%% Dimension reduction
for cur_nComponents = vNcomponents
    curDataParams.nComponents = cur_nComponents;

    MeISna = false;
    MlooISna = false;
    %% Majority Rule
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for w_n = 0:2
        curDataParams.wiggle_num = w_n;

        [mFileNam,~] = matFileNamSelect('pClassification',curDataParams);
        load([mFileNam matExt],'class_out','class_out_LOO')

        % Check for NA's
        if sum(isnan(class_out)) + sum(isnan(class_out_LOO)) > 0
            if sum(isnan(class_out)) > 0
                MeISna = true;
            end
            if sum(isnan(class_out_LOO)) > 0
                MlooISna = true;
            end
            break
        end

        if w_n == 0
            class_mr = class_out;
            class_mr_LOO = class_out_LOO;
        else
            class_mr = class_mr + class_out;
            class_mr_LOO = class_mr_LOO + class_out_LOO;
        end
    end
    
    class_mr = round(class_mr/curDataParams.wiggle_den);
    class_mr_LOO = round(class_mr_LOO/curDataParams.wiggle_den);
    
    curDataParams.dataType = 'Binary';
    curDataParams.normalisation = false;
    [mFileNam,~] = matFileNamSelect('pSummary',curDataParams);
    load([mFileNam matExt],'p_lnm','p_suit')    
    curDataParams.dataType = vDatTypes{datType_idx};
    curDataParams.normalisation = n;

    if MeISna
        Me_str = 'NA';
    else
        switch curDataParams.classificationMethod
            case {'NB' 'DWD'}
                Me_str = num2str(sum(class_mr(p_suit) ~= p_lnm(p_suit)'));
            otherwise
                Me_str = num2str(sum(class_mr ~= p_lnm(p_suit)'));                
        end
    end
    if MlooISna
        Mloo_str = 'NA';
    else
        Mloo_str = num2str(sum(class_mr_LOO ~= p_lnm(p_suit)'));
    end        

    fprintf(fid,[Me_str ',' ...
          Mloo_str ',' ...
          num2str(sum(p_suit)) ',' ...
          curDataParams.folNam ',' ...
          vr_method_cur ',' ...
          da_method_cur ',' ...
          curDataParams.dataType ',' ...
          num2str(curDataParams.normalisation) ',' ...
          num2str(curDataParams.includeEmptyVals) ',' ...
          num2str(curDataParams.smoothParam) ',' ...
          num2str(curDataParams.use_cancer_annot) ',' ...
          num2str(curDataParams.nCal_tol) ',' ...
          num2str(curDataParams.nSpec_r_tol) ',' ...
          num2str(curDataParams.nSpec_bin_tol) ',' ...
          num2str(curDataParams.nComponents) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
end    
end     
end
end
end
end
end
end
end
end
end
end

fclose(fid);













































