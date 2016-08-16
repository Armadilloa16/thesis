clc
clear all
close all

%%
% Add my code folders.
addpath('./functionMfiles/dataManagement',...
    './functionMfiles/analysisFunctions',...
    './functionMfiles/plottingFunctions')
% Some of Steve Marrons code, most particularly pcaSM for PCA.
addpath('./functionMfiles/codeSM/General',...
    './functionMfiles/codeSM/BatchAdjust',...
    './functionMfiles/codeSM/Smoothing')
% Some other peoples code I have used.
addpath('./functionMfiles/otherPeoplesCode')        
% addpath('./SDPT3-4.0/Solver')        

%%
[isOctave,matType,matExt] = checkIsOctave();
  
vNcal = [0];
% vNcal = [0 3 4];
v_minNspec_r = [1];
% v_minNspec_r = [1 100];

%%
% Calibrant Details.
% calibrantMZ = [1296.685,1570.677,2147.199,2932.588];
% calibrantName = {'Angiotensin I','[Glu]-Fibrinopeptide B', ... 
%     'Dynorphin A','ACTH fragment (1-24)'};

%%
% Proccessing parameters
curDataParams = struct();
curDataParams.binSize = 0.25;
curDataParams.ppmNormTol = true;
curDataParams.normTol = 10;


disp('----------------------')
disp('  Summarise Patients  ')
disp('----------------------')


% Datasets
dataGroupNam = 'Vtma';
vvFolNams = {{'VA1','VA2','VA3',...
              'VB1','VB2','VB3'}};

for vFolNams_idx = 1:length(vvFolNams)
vFolNams = vvFolNams{vFolNams_idx};    
disp(' ')
if length(vFolNams) > 1
    disp(dataGroupNam)
else
    disp(vFolNams{1})
end
    
%% Data Type
vDatTypes = {'Binary','LogIntensity','Intensity','SN','Area'};
for datType_idx = 1:length(vDatTypes) 
curDataParams.dataType = vDatTypes{datType_idx};
disp(['  ' curDataParams.dataType ' Data'])
switch curDataParams.dataType
    case 'Binary'
        curDataParams.normalisation = false;
        %% Binary Smoothing
        for curSmoothParam = [0 0.15 0.25]
        curDataParams.smoothParam = curSmoothParam;
        disp(['    ' num2str(curDataParams.smoothParam) '-smooth'])

        %% Wiggled Binning
        curDataParams.wiggle = true;
        curDataParams.wiggle_den = 3;
        for w_n = 0:(curDataParams.wiggle_den-1)
        curDataParams.wiggle_num = w_n;
        disp(['      Wiggle ' num2str(w_n) ' / ' num2str(curDataParams.wiggle_den)])    
        
        %% Annotated Cancer Spectra
        for use_cancer_annot = [false true]
        curDataParams.use_cancer_annot = use_cancer_annot;
        if use_cancer_annot
            disp(['        ' 'Using only annotated cancer spectra'])
        else
            disp(['        ' 'Using all spectra'])
        end

        %% Calibrant QA filter
        for nCal_tol = vNcal
        curDataParams.nCal_tol = nCal_tol;
        disp(['          ' num2str(curDataParams.nCal_tol) ' calibrant QA filter'])

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            vC = cell(length(vFolNams),1);
            p_number_temp = [];
            r_list_count_temp = [];
            dataset_idx = [];
            for folnam_idx = 1:length(vFolNams)
            curDataParams.folNam = vFolNams{folnam_idx};

                % Initial review of datasets
                [mFileNam,~] = matFileNamSelect('rSummary',curDataParams);
                load([mFileNam matExt],'r_list','r_list_count')
                r_list_count_temp = [r_list_count_temp r_list_count];
                dataset_idx = [dataset_idx repmat(folnam_idx,1,length(r_list))];
                [mFileNam,~] = matFileNamSelect('Annotation',curDataParams);
                load([mFileNam matExt],'p_number')
                p_number = p_number(r_list);
                p_number_temp = [p_number_temp p_number];
                [mFileNam,~] = matFileNamSelect('Binned',curDataParams);
                load([mFileNam matExt],'vbincentrs')
                vC{folnam_idx} = vbincentrs;

            end
            p_number = p_number_temp;
            r_list_count = r_list_count_temp;
            vC = combineDims(vC);
            clear *temp
                    %%%%% p_number        : Patient #'s for combined data
                    %%%%% dataset_idx     : Reversability information
                    %%%%% r_list_count    : # spectra per region
                    %%%%% vC              : Variable matching between datasets.

            nSpec_temp = zeros(size(vC,2),1);
            rdata_temp = zeros(size(vC,2),length(p_number));
            for folnam_idx = 1:length(vFolNams)
            curDataParams.folNam = vFolNams{folnam_idx};

                [mFileNam,~] = matFileNamSelect('rSummary',curDataParams);
                load([mFileNam matExt],'rdata','nSpec')
                rdata_temp(vC(folnam_idx,:)' > 0,dataset_idx == folnam_idx) = rdata;
                nSpec_temp(vC(folnam_idx,:)' > 0) = nSpec_temp(vC(folnam_idx,:)' > 0) + nSpec;

            end
            rdata = rdata_temp;
            nSpec = nSpec_temp;
            clear *temp
                    %%%%% rdata        : Combined Region Summarised Data
                    %%%%% nSpec        : Combined # spectra with peaks in each bin (total)

                    
            for minNspec_r = v_minNspec_r
            curDataParams.nSpec_r_tol = minNspec_r;
            disp(['            ' 'minNspecPerRegion: ' num2str(minNspec_r)])
            r_subset = r_list_count >= curDataParams.nSpec_r_tol;
                    %%%%% r_subset     : Regions that satisfy min # spectra filter.

            % Load Clinical Data        
            [mFileNam,~] = matFileNamSelect('Vtma Clinical Data');
            load([mFileNam matExt],'p_list','p_lnm','p_suit')
            temp = zeros(max(p_list),1);
            temp(p_list) = p_lnm;
            p_lnm = temp == 1;
            temp(p_list) = p_suit;
            p_suit = temp == 1;
            temp(p_list) = p_list;
            p_list = temp;
            clear temp

            p_list_temp = unique(p_number(p_number~=0 & r_subset));
            p_lnm  = p_lnm(p_list_temp);
            p_suit = p_suit(p_list_temp);
            p_list = p_list(p_list_temp);
            clear p_list_temp
                    %%%%% p_list     : List of patients remaining after preprocessing.
                    %%%%% p_lnm      : LNM status of remaining patients
                    %%%%% p_suit     : Suitability of remaining patients

            nVar = size(rdata,1);
            pdata = zeros(nVar,length(p_list));
            for p = p_list'
                s = r_subset & p_number == p;
                pdata(:,p_list==p) = sum(rdata(:,s).*repmat(r_list_count(s),nVar,1),2)/sum(r_list_count(s));
            end

            % Save Results
            if length(vFolNams) > 1
                curDataParams.folNam = dataGroupNam;
            end
            vbincentrs = vC(end,:)';
            [mFileNam,regexpVars] = matFileNamSelect('pSummary',curDataParams);
            save(matType,[mFileNam matExt],'-regexp',regexpVars)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        end
        end
        end
        end
        
    case {'LogIntensity','Intensity','SN','Area'}
        
        curDataParams.smoothParam = 0;
        %% Normalisation
        for n = [false true]
        curDataParams.normalisation = n;
        if n
            disp('    With normalisation')
        else
            disp('    No normalisation')
        end

        %% Wiggled Binning
        curDataParams.wiggle = true;
        curDataParams.wiggle_den = 3;
        for w_n = 0:(curDataParams.wiggle_den-1)
        curDataParams.wiggle_num = w_n;
        disp(['      Wiggle ' num2str(w_n) ' / ' num2str(curDataParams.wiggle_den)])    

        %% Annotated Cancer Spectra
        for use_cancer_annot = [false true]
        curDataParams.use_cancer_annot = use_cancer_annot;
        if use_cancer_annot
            disp(['        ' 'Using only annotated cancer spectra'])
        else
            disp(['        ' 'Using all spectra'])
        end
        
        %% Calibrant QA filter
        for nCal_tol = vNcal
        curDataParams.nCal_tol = nCal_tol;
        disp(['          ' num2str(curDataParams.nCal_tol) ' calibrant QA filter'])
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            vC = cell(length(vFolNams),1);
            p_number_temp = [];
            r_list_count_temp = [];
            dataset_idx = [];
            curDataParams.dataType = 'Binary';
            curDataParams.normalisation = false;
        
            for folnam_idx = 1:length(vFolNams)
            curDataParams.folNam = vFolNams{folnam_idx};

                % Initial review of datasets
                [mFileNam,~] = matFileNamSelect('rSummary',curDataParams);
                load([mFileNam matExt],'r_list','r_list_count')
                r_list_count_temp = [r_list_count_temp r_list_count];
                dataset_idx = [dataset_idx repmat(folnam_idx,1,length(r_list))];
                [mFileNam,~] = matFileNamSelect('Annotation',curDataParams);
                load([mFileNam matExt],'p_number')
                p_number = p_number(r_list);
                p_number_temp = [p_number_temp p_number];
                [mFileNam,~] = matFileNamSelect('Binned',curDataParams);
                load([mFileNam matExt],'vbincentrs')
                vC{folnam_idx} = vbincentrs;
                                
            end
            p_number = p_number_temp;
            r_list_count = r_list_count_temp;
            vC = combineDims(vC);
            curDataParams.dataType = vDatTypes{datType_idx};
            curDataParams.normalisation = n;
            clear *temp
                    %%%%% p_number        : Patient #'s for combined data
                    %%%%% dataset_idx     : Reversability information
                    %%%%% r_list_count    : # spectra per region
                    %%%%% vC              : Variable matching between datasets.

            nSpec_temp = zeros(size(vC,2),1);
            rdata_temp = zeros(size(vC,2),length(p_number));
            for folnam_idx = 1:length(vFolNams)
            curDataParams.folNam = vFolNams{folnam_idx};

                curDataParams.dataType = 'Binary';
                curDataParams.normalisation = false;
                [mFileNam,~] = matFileNamSelect('rSummary',curDataParams);
                load([mFileNam matExt],'nSpec')
                curDataParams.dataType = vDatTypes{datType_idx};
                curDataParams.normalisation = n;
                [mFileNam,~] = matFileNamSelect('rSummary',curDataParams);
                load([mFileNam matExt],'rdata')
                rdata_temp(vC(folnam_idx,:)' > 0,dataset_idx == folnam_idx) = rdata;
                nSpec_temp(vC(folnam_idx,:)' > 0) = nSpec_temp(vC(folnam_idx,:)' > 0) + nSpec;

            end
            rdata = rdata_temp;
            nSpec = nSpec_temp;
            clear *temp
                    %%%%% rdata        : Combined Region Summarised Data
                    %%%%% nSpec        : Combined # spectra with peaks in each bin (total)

            for minNspec_r = v_minNspec_r
            curDataParams.nSpec_r_tol = minNspec_r;
            disp(['            ' 'minNspecPerRegion: ' num2str(minNspec_r)])
            r_subset = r_list_count >= curDataParams.nSpec_r_tol;
                    %%%%% r_subset     : Regions that satisfy min # spectra filter.

            % Load Clinical Data        
            [mFileNam,~] = matFileNamSelect('Vtma Clinical Data');
            load([mFileNam matExt],'p_list','p_lnm','p_suit')
            temp = zeros(max(p_list),1);
            temp(p_list) = p_lnm;
            p_lnm = temp == 1;
            temp(p_list) = p_suit;
            p_suit = temp == 1;
            temp(p_list) = p_list;
            p_list = temp;
            clear temp

            p_list_temp = unique(p_number(p_number~=0 & r_subset));
            p_lnm  = p_lnm(p_list_temp);
            p_suit = p_suit(p_list_temp);
            p_list = p_list(p_list_temp);
            clear p_list_temp
                    %%%%% p_list     : List of patients remaining after preprocessing.
                    %%%%% p_lnm      : LNM status of remaining patients
                    %%%%% p_suit     : Suitability of remaining patients

            nVar = size(rdata,1);
            pdata = zeros(nVar,length(p_list));
            for p = p_list'
                s = r_subset & p_number == p;
                pdata(:,p_list==p) = sum(rdata(:,s).*repmat(r_list_count(s),nVar,1),2)/sum(r_list_count(s));
            end

            % Save Results
            if length(vFolNams) > 1
                curDataParams.folNam = dataGroupNam;
            end
            vbincentrs = vC(end,:)';
            [mFileNam,regexpVars] = matFileNamSelect('pSummary',curDataParams);
            save(matType,[mFileNam matExt],'-regexp',regexpVars)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        end
        end
        end
        end
        
        
end
end
end



        
        





