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

vUseAnnot = [false true];

vNcal = 0;
% vNcal = [0 3 4];
v_minNspec_r = 1;
% v_minNspec_r = [1 100];
v_minNspec_bin = 1;
% v_minNspec_bin = [1 100];
vCCArankTypes = 2;
% vCCArankTypes = 1:2;

fid = fopen('.\output\Vtma_cca_ranked_variables.txt','wt');

v_cols = 1:30;

% Write Header
fprintf(fid,['VRmethod,DataType,Normalisation,includeEmptyValues,' ...
    'Smooth,CancerAnnotation,minNcal,minNspecPerCore,minNspecPerBin,' ...
    'Wiggle' num2str(v_cols(1:9),',V%i') ...
    num2str(v_cols(10:end),',V%i') '\n']);

%% Dataset
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};

%% Data Type
for datType_idx = 1:length(vDatTypes) 
curDataParams.dataType = vDatTypes{datType_idx};
disp([curDataParams.dataType ' Data'])

%% Annotated Cancer Spectra
for use_cancer_annot = [false true]
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

%% Bin Wiggling
for w_n = 0:2
    curDataParams.wiggle_num = w_n;
    disp(['  ' '  ' num2str(curDataParams.wiggle_num) ' / ' num2str(curDataParams.wiggle_den) ' wiggle'])

    %% Load Patient Summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~strcmp(curDataParams.dataType,'Binary')
        curDataParams.dataType = 'Binary';
        curDataParams.normalisation = false;

        [mFileNam,regexpVars] = matFileNamSelect('pSummary',curDataParams);
        load([mFileNam matExt],'-regexp',regexpVars)
        pbdata = pdata;
        % Number of spectra per bin filter
        pbdata = pbdata(nSpec >= curDataParams.nSpec_bin_tol,:);

        curDataParams.dataType = vDatTypes{datType_idx};
        curDataParams.normalisation = n;
    end

    [mFileNam,regexpVars] = matFileNamSelect('pSummary',curDataParams);
    load([mFileNam matExt],'-regexp',regexpVars)

    % Number of spectra per bin filter
    pdata = pdata(nSpec >= curDataParams.nSpec_bin_tol,:);

    if ~curDataParams.includeEmptyVals
        pdata(pbdata > 0) = pdata(pbdata > 0)./pbdata(pbdata > 0);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    curDataParams.restrict_p_suit = true;
    [mFileNam,~] = matFileNamSelect('pPCA',curDataParams);
    load([mFileNam matExt],'mpc','meigvec','veigval')

    X = pdata(:,p_suit);
    Y = double(p_lnm(p_suit))';
    Z = meigvec*diag(veigval.^(-1/2))*meigvec';
    % Note that p1 is proportional to C in this case (because Y is one-dimensional)
    p1 = Z*(1/sum(p_suit))*(X-repmat(mean(X,2),1,size(X,2)))*(Y-mean(Y))'*(var(Y)^(-1/2));
    phi1 = Z*p1;

    for rank_type = vCCArankTypes
        switch rank_type
            case 1
                curDataParams.classificationMethod = 'altCCA';
                [~,I] = sort(abs(p1),'descend');
            case 2
                curDataParams.classificationMethod = 'CCA';
                [~,I] = sort(abs(phi1),'descend');
        end    
        
        v_ranks = num2str(vbincentrs(I(1:30)),',%4.4f')';
        v_ranks = v_ranks(:)';
        v_ranks = v_ranks(v_ranks ~= ' ');
        
        fprintf(fid,[curDataParams.classificationMethod ',' ...
            curDataParams.dataType ',' ...
            num2str(curDataParams.normalisation) ',' ...
            num2str(curDataParams.includeEmptyVals) ',' ...
            num2str(curDataParams.smoothParam) ',' ...
            num2str(curDataParams.use_cancer_annot) ',' ...
            num2str(curDataParams.nCal_tol) ',' ...
            num2str(curDataParams.nSpec_r_tol) ',' ...
            num2str(curDataParams.nSpec_bin_tol) ',' ...
            num2str(curDataParams.wiggle_num) ...
            v_ranks '\n']);
        
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













































