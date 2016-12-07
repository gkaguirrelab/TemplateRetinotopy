%% Setup the directories and variables
savedir = ['/Users/abock/Dropbox (Aguirre-Brainard Lab)/Bock/Bock_Aguirre_manuscripts/' ...
    'Retinotopy_Templates/Figures/Figure5_Partition/raw'];
cd(savedir);

sessions = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014/' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014'};
allRuns = {'[3,4,6]' '[2,4,6]' '[2,4,6]'}; % must be a string (!)
%templateType = 'coarse';
templateType = 'coarseV2V3size';
tcPart = 'half';
func = 's5.wdrf.tf';
V2V3 = '1';

logDir = '/data/jet/abock/LOGS';
hemis = {'lh' 'rh'};

fitType = 'V2V3';
badAreas = [4,5]; % visual areas to exclude;
eLims = [1 5]; % eccentricity lims [low high];

tDirName = 'coarse_model_templates_V2V3size';

%% Create the regress template scripts, for 1st and 2nd Halves
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    script_dir = fullfile(session_dir,'fit_template_scripts','partition');
    system(['rm -rf ' script_dir]);
    mkdir(script_dir);
    runs = allRuns{ss};
    for i = [0,3,5] % 30, 15, and 5 minute partitions
        tmp =combnk(1:6,i);
        if isempty(tmp)
            leaveOut = 0;
            leaveOutName = '';
            outDir = fullfile(script_dir,[templateType '_' tcPart '_leaveOut' leaveOutName]);
            saveDir = fullfile(session_dir,'pRFs',templateType,func,...
                ['Movie_' tcPart '_leaveOut' leaveOutName]);
            create_regress_template_scripts(session_dir,templateType,outDir,runs,...
                func,saveDir,tcPart,['[' num2str(leaveOut) ']'],V2V3);
        else
            for j = 1:size(tmp,1)
                leaveOut = tmp(j,:);
                leaveOutName = strrep(num2str(leaveOut),' ','');
                outDir = fullfile(script_dir,[templateType '_' tcPart '_leaveOut' leaveOutName]);
                saveDir = fullfile(session_dir,'pRFs',templateType,func,...
                    ['Movie_' tcPart '_leaveOut' leaveOutName]);
                create_regress_template_scripts(session_dir,templateType,outDir,runs,...
                    func,saveDir,tcPart,['[' num2str(leaveOut) ']'],V2V3);
            end
        end
    end
end
%% Submit the 5-min partition regress jobs
% Run this from chead
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    script_dir = fullfile(session_dir,'fit_template_scripts','partition');
    cDirs = listdir(fullfile(script_dir,[templateType '_' tcPart '_leaveOut*']),'dirs');
    for i = 1:length(cDirs);
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            system(['rm ' fullfile(logDir,'*')]);
            pause(5);
            scriptDir = fullfile(script_dir,cDirs{i});
            cd(scriptDir);
            inScript=['submit_' hemi '_regress.sh'];
            system(['sh ' fullfile(scriptDir,inScript)]);
            disp(['running ' fullfile(scriptDir,inScript)]);
            pause(5);
            system(['qstat > ' fullfile(logDir,'tmp.txt')]);
            pause(5);
            fid = fopen(fullfile(logDir,'tmp.txt'));
            tmp = fread(fid);
            fclose(fid);
            while ~isempty(tmp)
                system(['qstat > ' fullfile(logDir,'tmp.txt')]);
                pause(5);
                fid = fopen(fullfile(logDir,'tmp.txt'));
                tmp = fread(fid);
                fclose(fid);
            end
            fclose('all');
            rerun_regress(inScript,scriptDir,logDir);
            pause(5);
            system(['qstat > ' fullfile(logDir,'tmp.txt')]);
            pause(5);
            fid = fopen(fullfile(logDir,'tmp.txt'));
            tmp = fread(fid);
            fclose(fid);
            while ~isempty(tmp)
                system(['qstat > ' fullfile(logDir,'tmp.txt')]);
                pause(5);
                fid = fopen(fullfile(logDir,'tmp.txt'));
                tmp = fread(fid);
                fclose(fid);
            end
            fclose('all');
        end
    end
end
%% Compute partition error
splitComb = combnk(1:6,3);
partError = nan(length(sessions),63,length(hemis));
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    modelTdir = fullfile(session_dir,'pRFs','coarse_model_templates');
    aDirs = fullfile(session_dir,'pRFs',templateType,func);
    pDirs = listdir(fullfile(aDirs,'Movie_half*'),'dirs');
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        progBar = ProgressBar(length(pDirs),'computing error...');
        for p = 1:length(pDirs)
            tDir = fullfile(aDirs,pDirs{p});
            [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,[],[],[],fitType);
            dotinds = strfind(sorted_templates{1},'.');
            tempEcc = fullfile(modelTdir,...
                [hemi '.ecc.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']);
            tempPol = fullfile(modelTdir,...
                [hemi '.pol.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']);
            clear tmp
            for i = 1:length(splitComb)
                inEcc = fullfile(session_dir,'pRFs',...
                    [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                inPol = fullfile(session_dir,'pRFs',...
                    [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
                tmp(i) = nanmean(degerror);
            end
            partError(ss,p,hh) = mean(tmp(:));
            progBar(p);
        end
    end
end
save(fullfile(savedir,'partError.mat'),'partError');
%% Get the number of partitions in each partition directory
session_dir = sessions{1}; % It's the same for all session_dir
aDirs = fullfile(session_dir,'pRFs',templateType,func);
pDirs = listdir(fullfile(aDirs,'Movie_half*'),'dirs');
for p = 1:length(pDirs)
    tdir = fullfile(aDirs,pDirs{p});
    lInd = strfind(tdir,'leaveOut');
    numParts(p) = 6 - (length(tdir) - (lInd+7)); % (lInd+7) = length(leaveOut);
end
%% Get Means
for i = 1:max(numParts) % Number of partitions
    ct = 0;
    for ss = 1:length(sessions)
        for hh = 1:length(hemis)
            ct = ct + 1;
            tmp = partError(ss,numParts==i,hh);
            pMeans(i,ct) = mean(tmp(:));
        end
    end
end
%% Pool across subjects/hemis
for i = 1:max(numParts)
    tmp = pMeans(i,:);
    sError.mean(i) = mean(tmp(:));
    sError.SD(i) = std(tmp(:));
end
%% Error
fullFigure;
errorbar(1:max(numParts),sError.mean,sError.SD);
%plot(allParts,allMeans,'.','MarkerSize',20);
%hold on;
%plot(binParts,binMeans);
ylim([0 7]);
xlabel('Number of Partitions','FontSize',20);
ylabel('Mean error (degrees visual angle)','FontSize',20);
set(gca,'XTick',[1 2 3 4 5 6]);
%set(gca,'YTick',[0.05 0.1 0.15 0.2 0.25]);
axis square;
legend('Error bars = SD','Location','NorthEast');
savefigs('pdf','Partitions_vs_Error');
close all;
%% Compute error in split-halves
splitComb = combnk(1:6,3);
ct = 0;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        clear tmp
        for i = 1:length(splitComb)/2
            inEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            inPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            tempEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(end-(i-1),:),'%1d') '.prfs.nii.gz']);
            tempPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(end-(i-1),:),'%1d') '.prfs.nii.gz']);
            [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
            tmp(i) = nanmean(degerror);
        end
        splitError.mean(ct) = mean(tmp(:));
    end
end
save(fullfile(savedir,'splitError.mat'),'splitError');
%% Compute error in anat templates (deformed)
splitComb = combnk(1:6,3);
ct = 0;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    modelTdir = fullfile(session_dir,'pRFs','anat_templates');
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        tempEcc = fullfile(modelTdir,...
            [hemi '.ecc.anat.nii.gz']);
        tempPol = fullfile(modelTdir,...
            [hemi '.pol.anat.nii.gz']);
        clear tmp
        for i = 1:length(splitComb)
            inEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            inPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
            tmp(i) = nanmean(degerror);
        end
        anatsplitError.mean(ct) = mean(tmp(:));
    end
end
save(fullfile(savedir,'anatsplitError.mat'),'anatsplitError');
%% Compute error in 'best' coarse templates
splitComb = combnk(1:6,3);
ct = 0;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    tDir = fullfile(session_dir,'pRFs',templateType,func,'Movie',fitType);
    modelTdir = fullfile(session_dir,'pRFs',tDirName);
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,fitType);
        dotinds = strfind(sorted_templates{1},'.');
        tempEcc = fullfile(modelTdir,...
            [hemi '.ecc.' sorted_templates{1}(dotinds(1)+1:dotinds(6)-1) '.nii.gz']);
        tempPol = fullfile(modelTdir,...
            [hemi '.pol.' sorted_templates{1}(dotinds(1)+1:dotinds(6)-1) '.nii.gz']);
        clear tmp
        for i = 1:length(splitComb)
            inEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            inPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
            tmp(i) = nanmean(degerror);
        end
        coarsesplitError.mean(ct) = mean(tmp(:));
    end
end
%save(fullfile(savedir,'coarsesplitError.mat'),'coarsesplitError');
%% Compute error in anat templates (non-deformed)
splitComb = combnk(1:6,3);
ct = 0;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    modelTdir = fullfile(session_dir,'pRFs',tDirName);
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        tempEcc = fullfile(modelTdir,...
            [hemi '.ecc.4.4.4.nii.gz']);
        tempPol = fullfile(modelTdir,...
            [hemi '.pol.4.4.4.nii.gz']);
        clear tmp
        for i = 1:length(splitComb)
            inEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            inPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
            tmp(i) = nanmean(degerror);
        end
        nondanatsplitError.mean(ct) = mean(tmp(:));
    end
end
save(fullfile(savedir,'nondanatsplitError.mat'),'nondanatsplitError');
%% Plot error
dim = [.5 .75 .1 .1];
split.mean = mean(splitError.mean);
split.std = std(splitError.mean);
anat.mean = mean(anatsplitError.mean);
anat.std = std(anatsplitError.mean);
coarse.mean = mean(coarsesplitError.mean);
coarse.std = std(coarsesplitError.mean);
nodanat.mean = mean(nondanatsplitError.mean);
nodanat.std = std(nondanatsplitError.mean);
fullFigure;bar([split.mean anat.mean coarse.mean nodanat.mean]);
hold on;
errorbar([split.mean anat.mean coarse.mean nodanat.mean],...
    [split.std anat.std coarse.std nodanat.std],'.','MarkerSize',0.01);
xlabel('Map Type','FontSize',20);
ylabel('Mean error (degress visual angle)','FontSize',20);
set(gca,'XTickLabel',{'split-half','anat-deformed','coarse','anat-non-deformed'},'FontSize',15);
axis square
annotation('textbox',dim,'String','error bars = SD','FontSize',15,'HorizontalAlignment','center');
title('Error in retinotopy prediction','FontSize',20);
savefigs('pdf','templates_error');
close all;