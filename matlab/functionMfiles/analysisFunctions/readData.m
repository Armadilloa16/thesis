function [L,LXY,XYL,X,Y,fExists,Vars,emptySpec,R] = readData(folnam)
%%% INPUT PARAMETERS
%%% folnam:     folder containing the peaklists
%%% %%%
%%% OUTPUT PARAMETERS
%%% L:          Peak Lists as a nPeak_ul x nVars x nSpec_ul 3-array where
        %%% - nSpec_ul: The number of spectra
        %%% - nPeaks_ul: The maximum number of peaks in any given spectrum
        %%% - nVars: The number of things measured on each peak
%%% LXY:        A nSpec_ul x 3 array where
        %%% - The first column contains the X-index of the
        %%% corresponding spectrum
        %%% - The second column contains the Y-index of the corresponding
        %%% spectrum
        %%% - The third column contains the number of peaks featured in the
        %%% corresponding spectrum.
%%% XYL:        A 2D array index'd by X-index and Y-index, whose entries are zero
        %%% if no spectrum exists at that point, and by the index of that spectrum
        %%% in the list representations above if it does.
%%% X:          A vector of X-values, at their corresponding X-indices
%%% Y:          A vector of Y-values, at their corresponding Y-indices
%%% fExists:    A binary 2D array index'd by X-index and Y-index giving the
        %%% presence or absence of a spectrum at each point in space.
%%% Vars:       The Variable names of the things measured on each peak.

	more off

    errorLogColLabels = {'X' 'Y' 'm/z'};
    errorLog = zeros(0,3);
    errorType = {};
    
    dispProgress = true(1,1);
    Vars = {'m/z'	'SN'	'QualityFactor'	'Resolution'	'intensity'	'area'};
    nVars = length(Vars);

    %%% Find non-empty spectra
    if dispProgress
		disp('Compiling spatial map...')
    end
    % Check if peaklists exists?
%     fNames = ls(['./' folnam '/peaklists_quad']);
    fNames = ls(['./' folnam '/peaklists']);
%     fNames = ls(['./' folnam]);

    fNames = fNames';
    specNames = regexp(fNames(true(size(fNames)))','R\d{2}X\d{3}Y\d{3}','match');
    nSpec = length(specNames);
    if nSpec == size(fNames,2) - 2
        fNames = cellstr(fNames(:,3:end)');
    else
        error('There are non-peaklist files in the /peaklists folder.')
    end
    specNames = char(specNames)';
    Y_spec = [specNames(9:11,:); repmat(' ',1,nSpec)];
    Y_spec = str2num(Y_spec(true(size(Y_spec)))');
    X_spec = [specNames(5:7,:); repmat(' ',1,nSpec)];
    X_spec = str2num(X_spec(true(size(X_spec)))');
    R_spec = [specNames(2:3,:); repmat(' ',1,nSpec)];
    R_spec = str2num(R_spec(true(size(R_spec)))');
    clear specNames

    X = min(X_spec):max(X_spec);
    Y = min(Y_spec):max(Y_spec);
    fExists = false(length(X),length(Y));
    emptySpec = false(length(X),length(Y));
    XYL = zeros(size(fExists));
    for idx = 1:nSpec
        fExists(X_spec(idx) - X(1) + 1,Y_spec(idx) - Y(1) + 1) = true;
        XYL(X_spec(idx) - X(1) + 1,Y_spec(idx) - Y(1) + 1) = idx;
    end
    if sum(sum(fExists)) ~= nSpec
        error('nSpec does not match number of fExists')
    end
    LXY = zeros(nSpec,3);
    LXY(:,1) = X_spec(XYL(fExists))';
    LXY(:,2) = Y_spec(XYL(fExists))';
    R = R_spec(XYL(fExists))';
    
    % Correction for Region Numbers starting at zero.
    R = R + 1;

    fNames = fNames(XYL(fExists));
    
    XYL = zeros(size(fExists));
    XYL(fExists) = 1:sum(sum(fExists));
    L = cell(1,nSpec);
	if dispProgress
		disp('Reading in peaklists...')
    end
    progressCount = 1;
    for spec_idx = 1:nSpec
        
%         fname = [folnam '/peaklists_quad/' fNames{spec_idx}];
        fname = [folnam '/peaklists/' fNames{spec_idx}];
%         fname = [folnam '/' fNames{spec_idx}];
        
        [nPeaks,tempData] = readflat(fname,false);                
        sizeTempData = size(tempData);
        LXY(spec_idx,3) = nPeaks;
        
        if nPeaks == 0 % empty
            emptySpec(LXY(spec_idx,1) - X(1) + 1,LXY(spec_idx,2) - Y(1) + 1) = true;
            fExists(LXY(spec_idx,1) - X(1) + 1,LXY(spec_idx,2) - Y(1) + 1) = false;
        elseif sizeTempData(2) > nVars
            error(['ERROR: !!! number of variables is more than expected X:' num2str(LXY(spec_idx,1)) ' Y:' num2str(LXY(spec_idx,2))])
        elseif sizeTempData(2) < nVars
            error(['ERROR: !!! number of variables is less than expected X:' num2str(LXY(spec_idx,1)) ' Y:' num2str(LXY(spec_idx,2))])
        else
            if sum(sum(tempData == 0)) > 0
                disp(['WARNING: array contains zero values X:' num2str(LXY(spec_idx,1)) ' Y:' num2str(LXY(spec_idx,2))])
                tempI = find(sum(tempData == 0) > 0);
                tempDisp = '';
                for i = tempI
                    tempDisp = [tempDisp Vars{i} ' (' num2str(sum(tempData(:,i)==0)) ')   '];
                    
                    errorType = [errorType ; repmat({[Vars{i} ' zero']},sum(tempData(:,i) == 0),1)];
                    errorLog = [errorLog ; [repmat(LXY(spec_idx,1),sum(tempData(:,i) == 0),1) repmat(LXY(spec_idx,2),sum(tempData(:,i) == 0),1) tempData((tempData(:,i) == 0),1)]];
                    
                end
                disp(['     ' tempDisp])
                clear tempDisp tempI
            end
            L{spec_idx} = tempData;
			if dispProgress
                if spec_idx/nSpec >= progressCount/10;
					disp([num2str(progressCount*10) '% complete.'])
					progressCount = progressCount + 1;
                end
            end
        end
    end
    nonEmptySpec = XYL(fExists);
    XYL = zeros(size(fExists));
    XYL(fExists) = 1:sum(sum(fExists));
    L = L(nonEmptySpec);
    LXY = LXY(nonEmptySpec,:);
    
    
    
    save(['./errorLogMATfiles/errorLog_' folnam],'errorLogColLabels','errorLog','errorType')
                    
    
    
    
    
    %%%
    % Old shitty code
    %%% 
    
    
%     Xill = 0;
%     Xiul = 1000;
%     Yill = 0;
%     Yiul = 1000;
%     X = Xill:Xiul;
%     Y = Yill:Yiul;
% 
%     fExists = false(length(X),length(Y));
%     emptySpec = false(length(X),length(Y));
% 
%     nSpec_iul = 40000;
% 
%     LXY = zeros(nSpec_iul,3);
% 
%     count = 0;
% 	for j = 1:length(Y)
%         for i = 1:length(X)
%             fname = fnamSelect(folnam,X(i),Y(j));
%             fid = fopen(fname);        
%             if fid == -1 % No File at this X-Y coordinate
%                 fExists(i,j) = false;
%             else
%                 
%                 nPeaks = readflat(fname,false);            
%                 if nPeaks > 0 % File is non-empty
%                     fExists(i,j) = true;
%                     fclose(fid);                        
%                     count = count + 1;
%                     LXY(count,3) = nPeaks;
%                 elseif nPeaks == 0 % Spectrum exhibits no peaks.
%                     fExists(i,j) = false;
%                     emptySpec(i,j) = true;
%                     %if dispProgress
% 					%	disp(['Spectrum exhibiting no peaks at X: ' num2str(X(i)) ' Y: ' num2str(Y(j))])
%                     %end
% 					fclose(fid);
%                 end
%             end
%         end
% 		if dispProgress
% 			if mod(j,100) == 0
% 				disp([num2str(j/10) '% complete.'])
% 			end
% 		end
% 	end
% 	
% 	if dispProgress
% 		disp(' ')
% 	end
% 	
%     nSpec_ul = count;
%     LXY = LXY(1:nSpec_ul,:);
% 
%     Xll = find(sum(fExists,2),1,'first');
%     Xul = find(sum(fExists,2),1,'last');
%     Yll = find(sum(fExists,1),1,'first');
%     Yul = find(sum(fExists,1),1,'last');
% 
%     X = X(Xll:Xul);
%     Y = Y(Yll:Yul);
%     fExists = fExists(Xll:Xul,Yll:Yul);
%     emptySpec = emptySpec(Xll:Xul,Yll:Yul);













%     XYL = zeros(length(X),length(Y));
% 
%     nPeaks_ul = max(LXY(:,3));
%     L = zeros(nPeaks_ul,nVars,nSpec_ul);
% 
% 	if dispProgress
% 		disp('Reading in peaklists...')
% 	end
%     count = 0;
% 	progressCount = 1;
% 	nSpec = sum(sum(fExists));
%     for j = 1:length(Y)
%         for i = 1:length(X)
%             if fExists(i,j)
% 
%                 count = count + 1;
%                 LXY(count,1) = i;
%                 LXY(count,2) = j;
%                 XYL(i,j) = count;
% 
%                 fname = fnamSelect(folnam,X(i),Y(j));
%                 [nPeaks,tempData] = readflat(fname,false);
%                 sizeTempData = size(tempData);
% 
%                 if nPeaks ~= LXY(count,3)
%                     error(['ERROR: !!! number of peaks present do not match between detection run and read in run at X:' num2str(X(i)) ' Y:' num2str(Y(j))])
%                 elseif sizeTempData(2) > nVars
%                     error(['ERROR: !!! number of variables is more than expected X:' num2str(X(i)) ' Y:' num2str(Y(j))])
%                 elseif sizeTempData(2) < nVars
%                     error(['ERROR: !!! number of variables is less than expected X:' num2str(X(i)) ' Y:' num2str(Y(j))])
%                 else
%                     if sum(sum(tempData == 0)) > 0
%                         disp(['WARNING: array contains zero values X:' num2str(X(i)) ' Y:' num2str(Y(j))])
%                         disp(['     ' Vars(sum(tempData == 0) > 0)])
%                     end
%                     L(1:sizeTempData(1),:,count) = tempData;
% 					if dispProgress
% 						if count/nSpec >= progressCount/10;
% 							disp([num2str(progressCount*10) '% complete.'])
% 							progressCount = progressCount + 1;
% 						end
% 					end
% 				end
%             end
%         end
% 		%if dispProgress
% 		%	disp([num2str(j) ' columns complete']);
% 		%end
% 	end
% 	
% 	if dispProgress
% 		disp(' ')
% 	end
% 
%     % LL = zeros(1,(nVars+1),nPeaks_ul,nSpec_ul);
%     % for spec = 1:nSpec_ul
%     %     for k = 1:LXY(spec,3)
%     %         LL(1,1,k,spec) = k;
%     %         LL(1,2:(nVars+1),k,spec) = L(k,:,spec);
%     %     end
%     % end
% 
% 	more on
	
end 