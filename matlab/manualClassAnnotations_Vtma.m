clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement')

[isOctave,matType,matExt] = checkIsOctave();

disp('---------------------')
disp('     Annotation      ')
disp('---------------------')

% Clinical Data
raw_annot = csvread('..\data\Vtma_annotation\Vtma_clinical_data_numeric.csv',1,0);
p_list = raw_annot(:,1);
p_lnm  = raw_annot(:,2);
p_suit = true(length(p_list),1);
p_grade= raw_annot(:,3);

[mFileNam,regexpVars] = matFileNamSelect('Vtma Clinical Data');
save(matType,[mFileNam matExt],'-regexp',regexpVars)



% Preproccessing parameters
curDataParams = struct();

% Dataset
vFolNams = {'VA1','VA2','VA3',...
            'VB1','VB2','VB3'};

for folnam_idx = 1:length(vFolNams);
curDataParams.folNam = vFolNams{folnam_idx};
disp(curDataParams.folNam)

annot_cancer = csvread(['..\data\Vtma_annotation\' curDataParams.folNam '_annotation.csv'],1,0);
annot_cancer = annot_cancer(annot_cancer(:,1),2)' == 1;

mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'R','LXY')
R = unique(R);

% load([mFileNam matExt],'R','LXY','fExists','emptySpec')
% tmp = fExists + emptySpec == 1; 
% LXY = [LXY R(fExists(tmp))];
% R = unique(R);
% X = R;
% Y = R;
% i = 1;
% for r = R'
%     X(i) = mean(LXY(LXY(:,4)==r,1));
%     Y(i) = mean(LXY(LXY(:,4)==r,2));
%     i = i + 1;
% end

switch curDataParams.folNam

    % TMA 1B
    case {'VA1','VA2','VA3'}
        p_number = [ 10531513,12520576,13668910,5509498,7504891,12739050 ...
                     3518893,12512060,12739050,8509939,8513011,13534768 ...
                     9500915,10531513,8502283,13669790,11534968,13534768 ...
                     12520576,0,11531450,1524352,0,12527429 ...
                     0,13668910,0,12527429,8513011,11534968 ...
                     7514119,6510692,3518893,5515560,11543685,5515560 ...
                     6510692,0,0,0,1524352,7504891 ...
                     13539105,0,11543685,0,8502283,7514119 ...
                     13439500,12512060,8509939,11531450,13439500,0 ...
                     4517906,4517906,13539105,5509498,13668910,9500915];
        
    % TMA 2B
    case {'VB1','VB2','VB3'}
        p_number = [ 0,13534768,7504891,13534768,7514119,8513011 ...
                     11534968,3518893,12512060,13669790,0,8502283 ...
                     12512060,4517906,12520576,11534968,8502283,12527429 ...
                     5509498,3518893,12527429,8509939,11534968,0 ...
                     12739050,10531513,7514119,13644190,5515560,5515560 ...
                     11531450,8509939,11543685,0,1524352,13539105 ...
                     7504891,6510692,13439500,11531450,0,4517906 ...
                     0,6510691,9500915,0,10531513,5509498 ...
                     12739050,12520576,0,13539105,13439500,12543685 ...
                     13644190,0,13669790,9500915,0,8513011];

    otherwise
        error('Manual class annotations have not been added for this dataset.')
        
end

switch curDataParams.folNam
    case 'VA1'
        acq_order = [1:9 55:60 10:18 46:54 19:27 37:45 28:36];
        missing_cores = [4 5 7 9 11 15 25 26 48 57 58];
        ambiguous_regions = [];        
        
    case 'VA2'
        acq_order = [3:9 10:12 1 13:18 27:-1:19 28:36 37:44 2 45 54:-1:46 55:60];
        missing_cores = [4 5 9 11 15 26 36 48 58];
        ambiguous_regions = [35];        
        
    case 'VA3'
        acq_order = [1:9:55 56:-9:2 3:9:57 58:-9:4 5:9:59 60:-9:6 52:-9:7 53:-9:8 54:-9:9];
        missing_cores = [4 7 9 10 11 25 26 34 36 39 41 48 50 58];
        ambiguous_regions = [2];        
        
    case 'VB1'
        acq_order = [1:9 60:-1:55 10:18 54:-1:46 19:27 36 45 44 35 34 42 43 33 32 41 40 31 30 39 38 29 28 37];
        missing_cores = [1 5 14 24 37 39 41 55];
        ambiguous_regions = [];        
        
    case 'VB2'
        acq_order = [1:60];
        missing_cores = [1 5 14 37];
        ambiguous_regions = [];        
        
    case 'VB3'
        acq_order = [1:9 18:-1:10 55:60 46:54 27:-1:19 28:31 33 32 34:36 45:-1:37];
        missing_cores = [5 14 20 34 37 41 55];
        ambiguous_regions = [];        
end

acq_order = acq_order(~ismember(acq_order,missing_cores));
p_number(ambiguous_regions) = 0;
p_number = p_number(acq_order);
clear acq_order missing_cores ambiguous_regions
if strcmp(curDataParams.folNam,'VA3')
     p_number = [p_number zeros(1,33)];
end



% figure()
% if strcmp(curDataParams.folNam,'VA3')
%     plot(LXY(:,2),LXY(:,1),'.')
% else
%     plot(LXY(:,1),LXY(:,2),'.')
% end    
% title(curDataParams.folNam)
% switch curDataParams.folNam
%     case 'VA1'
%         set(gca,'YDir','Reverse')        
%     case 'VA2'
%         set(gca,'XDir','Reverse')                
%     case 'VA3'
%         set(gca,'XDir','Reverse','YDir','Reverse')                        
%     case 'VB1'
%         set(gca,'XDir','Reverse')                        
%     case 'VB2'
%         set(gca,'YDir','Reverse')                        
%     case 'VB3'
%         set(gca,'XDir','Reverse')                
% end
% hold on
% if strcmp(curDataParams.folNam,'VA3')
%     for i = 1:length(R)
%         text(Y(i),X(i),num2str(R(i)))
%     end
% else
%     for i = 1:length(R)
%         text(X(i),Y(i),num2str(R(i)))
%     end
% end    
% % if strcmp(curDataParams.folNam,'VA3')
% %     for i = 1:length(R)
% %         text(Y(i),X(i),num2str(p_number(i)))
% %     end
% % else
% %     for i = 1:length(R)
% %         text(X(i),Y(i),num2str(p_number(i)))
% %     end
% % end    

[mFileNam,regexpVars] = matFileNamSelect('Annotation',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)

dlmwrite(['./output/' curDataParams.folNam '_r_p_numbers.csv'],p_number')

end



