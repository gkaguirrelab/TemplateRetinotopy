%% TemplateRetinotopy MASTER script
%
%	This script will preprocess and analyze data for the TemplateRetinotopy project
%
%   MRklar - v1.0.0

%% Set session and subject name
logDir = '/data/jet/abock/LOGS';
dataDir = '/data/jet/abock/data';
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
templateDir = fullfile(dataDir,'Retinotopy_Templates','fsaverage_sym');
sessions = {...
    fullfile(dataDir,'Retinotopy_Templates','AEK','10012014') ...
    fullfile(dataDir,'Retinotopy_Templates','ASB','10272014') ...
    fullfile(dataDir,'Retinotopy_Templates','GKA','10152014') ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
jobNames = {'A100114K' 'A102714B' 'G101514A'};
outNames = {'AEK' 'ASB' 'GKA'};
pRFruns = {[1,2,5] [1,3,5] [1,3,5]};
pRFfuncs = {'wdrf.tf' 's5.wdrf.tf'};
srcROIs = {'cortex'};
movieRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
volFunc = 's5.wdrf.tf';% 'wdrf.tf';
hemis = {'lh' 'rh'};
pRFmem = 30; % memory for cluster
pRFmaps = {'co' 'coecc' 'copol' 'copeakt' 'cosig1' 'cosig2' 'cosig3' 'cosig4'};
tDirName = 'coarse_model_templates_V2V3size';
%% Run preprocessing
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    job_name = jobNames{ss};
    
    outDir = fullfile(session_dir,'preprocessing_scripts');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    numRuns = 6; % number of bold runs
    reconall = 0;
    slicetiming = 1; % correct slice timings
    refvol = 1; % motion correct to 1st TR
    filtType = 'detrend';
    lowHz = 0.01;
    highHz = 0.10;
    physio = 1;
    motion = 1;
    task = 0;
    localWM = 1;
    anat = 1;
    amem = 20;
    fmem = 50;
    create_preprocessing_scripts(session_dir,subject_name,outDir,logDir,...
        job_name,numRuns,reconall,slicetiming,refvol,filtType,lowHz,highHz,...
        physio,motion,task,localWM,anat,amem,fmem);
end
%% Run the above shell scripts in terminal
%   for example:
%   sh /data/jet/abock/data/Retinotopy_Template/AEK/10012014/preprocessing_scripts/submit_A100114K_all.sh

%% Temporally filter the rf.nii.gz
func = 'rf';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    boldDirs = find_bold(session_dir);
    for i = 1:length(boldDirs)
        runNum = i;
        temporal_filter(session_dir,runNum,func,filtType);
    end
end
%% Project retinotopic templates to subject space
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    project_template(session_dir,subject_name)
end
%% Decimate the anatomical template files
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    tDir = fullfile(session_dir,'anat_templates');
    decimate_templates(subject_name,tDir);
end
%% Create pRF prediction variables
% params.fieldSize         = 6.2116; % radius of stimuluated visual field in degrees visual angle
% params.TR                = 2; % TR is seconds
% params.peakHRF           = 3:.5:10; % time to HRF peak
% params.nFrames           = 300; % number of images
% params.cntrSig           = [1/4 1/4 20]; % min/lin_step/max sigma in degrees visual angle
% params.sigBorder         = 5; % Linear scale for visual field less than 5 degrees
% params.nLogSigs          = 6; % number of log sigmas
% params.srndSig           = 0; % scale of surround sigma
% params.ctrBetas          = 1; % center amplitude
% params.srndBetas         = 0; % surround amplitude
% params.gridScale         = 2; % scale the size of the stimulated visual field
% createPredpRF(imFile,paramsFile,outFile,params)
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    outDir = fullfile(session_dir,'pRF_prediction_scripts');
    system(['rm -rf ' outDir]);
    mkdir(outDir);
    imFiles = cell(1,length(pRFruns{ss}));
    paramFiles = cell(1,length(pRFruns{ss}));
    outFiles = cell(1,length(pRFruns{ss}));
    for i = 1:length(pRFruns{ss})
        imFiles{i} = fullfile(session_dir,'Stimuli',['run' num2str(pRFruns{ss}(i))],'bars_images.mat');
        paramFiles{i} = fullfile(session_dir,'Stimuli',['run' num2str(pRFruns{ss}(i))],'bars_params.mat');
        outFiles{i} = fullfile(session_dir,'Stimuli',['run' num2str(pRFruns{ss}(i))],'pRFpred.mat');
    end
    create_pRF_prediction_scripts(outDir,logDir,imFiles,paramFiles,outFiles,15);
end
%% Create pRF analysis scripts
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    outDir = fullfile(session_dir,'pRF_calc_scripts');
    system(['rm -rf ' outDir]);
    mkdir(outDir);
    boldDirs = find_bold(session_dir);
    ct = 0;
    predFiles = cell(1,length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    outFiles = cell(1,length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    inFiles = cell(1,length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    srcInds = cell(1,length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    for i = 1:length(pRFruns{ss})
        for hh = 1:length(hemis)
            for ff = 1:length(pRFfuncs)
                ct = ct + 1;
                outFiles{ct} = fullfile(session_dir,'pRFs',boldDirs{pRFruns{ss}(i)},...
                    [hemis{hh} '.' pRFfuncs{ff} '.pRFcalc.mat']);
                predFiles{ct} = fullfile(session_dir,'Stimuli',...
                    ['run' num2str(pRFruns{ss}(i))],'pRFpred.mat');
                inFiles{ct} = fullfile(session_dir,boldDirs{pRFruns{ss}(i)},...
                    [pRFfuncs{ff} '.surf.' hemis{hh} '.nii.gz']);
                tmp = load_nifti(inFiles{ct});
                srcInds{ct} = ['1:' num2str(size(tmp.vol,1))];
            end
        end
    end
    create_pRF_calc_scripts(outDir,logDir,outFiles,predFiles,inFiles,srcInds,pRFmem);
end
%% Make split-halves volumes
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    boldDirs = find_bold(session_dir);
    for i = 1:length(pRFruns{ss})
        for hh = 1:length(hemis)
            for ff = 1:length(pRFfuncs)
                inFile = fullfile(session_dir,boldDirs{pRFruns{ss}(i)},...
                    [pRFfuncs{ff} '.surf.' hemis{hh} '.nii.gz']);
                outFile1 = fullfile(session_dir,boldDirs{pRFruns{ss}(i)},...
                    [pRFfuncs{ff} '.surf.' hemis{hh} '.split1.nii.gz']);
                outFile2 = fullfile(session_dir,boldDirs{pRFruns{ss}(i)},...
                    [pRFfuncs{ff} '.surf.' hemis{hh} '.split2.nii.gz']);
                tmp = load_nifti(inFile);
                out1 = tmp;
                out2 = tmp;
                out1.vol = out1.vol(:,:,:,1:150);
                out2.vol = out2.vol(:,:,:,151:300);
                save_nifti(out1,outFile1);
                save_nifti(out2,outFile2);
            end
        end
    end
end
%% Create pRF split-halves analysis scripts
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    outDir = fullfile(session_dir,'pRF_split_calc_scripts');
    system(['rm -rf ' outDir]);
    mkdir(outDir);
    boldDirs = find_bold(session_dir);
    ct = 0;
    predFiles = cell(1,2*length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    outFiles = cell(1,2*length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    inFiles = cell(1,2*length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    srcInds = cell(1,2*length(pRFruns{ss})*length(hemis)*length(pRFfuncs));
    for i = 1:length(pRFruns{ss})
        for hh = 1:length(hemis)
            for ff = 1:length(pRFfuncs)
                for j = 1:2
                    ct = ct + 1;
                    outFiles{ct} = fullfile(session_dir,'pRFs',boldDirs{pRFruns{ss}(i)},...
                        [hemis{hh} '.' pRFfuncs{ff} '.split' num2str(j) '.pRFcalc.mat']);
                    predFiles{ct} = fullfile(session_dir,'Stimuli',...
                        ['run' num2str(pRFruns{ss}(i))],'pRFpred.mat');
                    inFiles{ct} = fullfile(session_dir,boldDirs{pRFruns{ss}(i)},...
                        [pRFfuncs{ff} '.surf.' hemis{hh} '.split' num2str(j) '.nii.gz']);
                    tmp = load_nifti(inFiles{ct});
                    srcInds{ct} = ['1:' num2str(size(tmp.vol,1))];
                end
            end
        end
    end
    create_pRF_calc_scripts(outDir,logDir,outFiles,predFiles,inFiles,srcInds,pRFmem);
end
%% Save the pRFs on the cortex/volume
srcROI = 'cortex';
progBar = ProgressBar(length(sessions),'pRFing...');
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    boldDirs = find_bold(session_dir);
    for i = 1:length(pRFruns{ss})
        outDir = fullfile(session_dir,'pRFs',boldDirs{pRFruns{ss}(i)});
        matFiles = listdir(fullfile(outDir,'*pRFcalc.mat'),'files');
        for j = 1:length(matFiles)
            matFile = fullfile(outDir,matFiles{j});
            % Set output name
            nameInd = strfind(matFiles{j},'.pRFcalc');
            outName = matFiles{j}(1:nameInd-1);
            hemi = matFiles{j}(1:2);
            templateFile = fullfile(session_dir,'anat_templates',[hemi '.ecc.anat.nii.gz']);
            tmp = load_nifti(templateFile);
            srcind = 1:length(tmp.vol);
            savepRF(matFile,outName,outDir,srcROI,srcind,templateFile);
        end
    end
    progBar(ss);
end

%% Average pRF
srcROI = 'cortex';
progBar = ProgressBar(length(sessions),'Averaging pRFs...');
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    boldDirs = find_bold(session_dir);
    for rr = 1:length(pRFfuncs)
        rootName = pRFfuncs{rr};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            for mm = 1:length(pRFmaps)
                pRFmap = pRFmaps{mm};
                for i = 1:length(pRFruns{ss})
                    inFiles{i} = fullfile(session_dir,'pRFs',boldDirs{pRFruns{ss}(i)},...
                        [hemi '.' rootName '.' srcROI '.' pRFmap '.prfs.nii.gz']);
                end
                outName = fullfile(session_dir,'pRFs',...
                    [hemi '.' rootName '.' srcROI '.' pRFmap '.avg.prfs.nii.gz']);
                average_pRF(inFiles,outName,srcROI,pRFmap);
            end
        end
    end
    progBar(ss);
end
%% Project pRF to fsaverage_sym space (all project to left hemisphere)
srcROI = 'cortex';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    boldDirs = find_bold(session_dir);
    for rr = 1:length(pRFfuncs)
        rootName = pRFfuncs{rr};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            for mm = 1:length(pRFmaps)
                pRFmap = pRFmaps{mm};
                sval = fullfile(session_dir,'pRFs',...
                    [hemi '.' rootName '.' srcROI '.' pRFmap '.avg.prfs.nii.gz']);
                tval = fullfile(session_dir,'pRFs',...
                    [hemi '.' rootName '.' srcROI '.' pRFmap '.avg.prfs.sym.nii.gz']);
                if strcmp(hemi,'lh')
                    mri_surf2surf(subject_name,'fsaverage_sym',sval,tval,hemi);
                else
                    mri_surf2surf([subject_name '/xhemi'],'fsaverage_sym',sval,tval,'lh');
                end
            end
        end
    end
    progBar(ss);
end
%% Average pRF split-halves
srcROI = 'cortex';
splitComb = combnk(1:6,3);
progBar = ProgressBar(length(sessions),'Averaging pRFs...');
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    boldDirs = find_bold(session_dir);
    for rr = 1:length(pRFfuncs)
        rootName = pRFfuncs{rr};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            for mm = 1:length(pRFmaps)
                pRFmap = pRFmaps{mm};
                ct = 0;
                clear tmpFiles
                for i = 1:length(pRFruns{ss})
                    theseFiles = listdir(fullfile(session_dir,'pRFs',boldDirs{pRFruns{ss}(i)},...
                        [hemi '.' rootName '.split*.' srcROI '.' pRFmap '.prfs.nii.gz']),'files');
                    for j = 1:length(theseFiles)
                        ct = ct + 1;
                        tmpFiles{ct} = fullfile(session_dir,'pRFs',boldDirs{pRFruns{ss}(i)},...
                            theseFiles{j});
                    end
                end
                for sc = 1:length(splitComb)
                    inFiles = tmpFiles(splitComb(sc,:));
                    outName = fullfile(session_dir,'pRFs',...
                        [hemi '.' rootName '.' srcROI '.' pRFmap '.avg.' ...
                        num2str(splitComb(sc,:),'%1d') '.prfs.nii.gz']);
                    average_pRF(inFiles,outName,srcROI,pRFmap);
                end
            end
        end
    end
    progBar(ss);
end
%% Prepare pRF template for fitting in Mathematica
for ff = 1:length(pRFfuncs);
    for ss = 1:length(sessions)
        session_dir = sessions{ss};
        subject_name = subjects{ss};
        outName = outNames{ss};
        prepare_pRF_Mathematica(session_dir,subject_name,outName,pRFfuncs{ff});
    end
end
%% Run Mathematica template fitting
%
%
%% Create pRF templates following Mathematica template fitting
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    create_pRF_template(session_dir);
    tDir = fullfile(session_dir,'pRFs','pRF_templates');
    decimate_templates(subject_name,tDir);
end
%% Copy an anatomical template, to be used later
mkdir(fullfile(templateDir,'anat_templates'));
system(['cp ~/data/2014-10-29.areas-template.nii.gz ' ...
    fullfile(templateDir,'anat_templates','lh.areas.anat.nii.gz')]);
%% Convert fsaverage_sym templates to nifti (coarse)
mem = 10;
outDir = fullfile(templateDir,'convertTemplateScripts');
tDir = fullfile(templateDir,'pRFs',tDirName);
create_convert_Mathematica_template_scripts(outDir,logDir,templateDir,tDir,mem);
% convert_Mathematica_templates(templateDir,tDir);

%%% submit the scripts above %%%

%% Project fsaverage_sym templates to subject space (coarse)
% This will create scripts to project the templates from the fsaverag_sym
% surface to the individual subject surfaces
tDir = fullfile(templateDir,'pRFs',tDirName);
tFiles = listdir(fullfile(tDir,'*.nii.gz'),'files');
hemis = {'lh' 'rh'};
srcsubject = 'fsaverage_sym';
mem = 5;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    submit_name = jobNames{ss};
    outDir = fullfile(session_dir,[tDirName '_scripts']);
    system(['rm -rf ' outDir]);
    mkdir(outDir);
    sDir = fullfile(session_dir,'pRFs',tDirName);
    system(['rm -rf ' sDir]);
    mkdir(sDir);
    job_string = cell(1,length(tFiles)*2);
    ct = 0;
    for tt = 1:length(tFiles);
        sval = fullfile(tDir,tFiles{tt});
        for hh = 1:length(hemis)
            ct = ct +1;
            hemi = hemis{hh};
            job_name = [hemi tFiles{tt}(3:end)];
            tval = fullfile(sDir,job_name);
            if strcmp(hemi,'lh')
                matlab_string = ['mri_surf2surf(''' srcsubject ''',''' subject_name ''',''' ...
                    sval ''',''' tval ''',''lh'');'];
            else
                matlab_string = ['mri_surf2surf(''' srcsubject ''',''' subject_name '/xhemi'',''' ...
                    sval ''',''' tval ''',''lh'');'];
            end
            create_job_shell(outDir,job_name,matlab_string);
            job_string{ct} = fullfile(outDir,[job_name '.sh']);
        end
    end
    create_submit_shell(outDir,logDir,submit_name,job_string,mem)
end
%% Run the above shell scripts
%   for example:
%   sh /data/jet/abock/data/Retinotopy_Templates/GKA/10152014/coarse_model_templates_V2V3size_scripts/G101514A.sh

%% Decimate the templates (coarse)
for ss = 1:length(sessions)
    params.outDir       = fullfile(sessions{ss},['decimate_' tDirName '_scripts']);
    system(['rm -rf ' params.outDir]);
    mkdir(params.outDir);
    params.logDir       = logDir;
    params.subjectName  = subjects{ss};
    params.templateDir  = fullfile(sessions{ss},'pRFs',tDirName);
    params.mem          = 5;
    create_decimate_template_scripts(params);
end

%%% submit the scripts above %%%

%% Decimate the bold runs
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    decimate_bold(session_dir,subject_name,volFunc);
end
%% Create cluster shell scripts (pRF,anat,coarse)
%tTypes = {'anat' 'pRF' 'coarse' 'coarseV2V3size'};
tTypes = {'coarseV2V3size'};
fitTypes = {'V1' 'V2V3'};
for ff = 1:length(fitTypes)
    for tt = 1:length(tTypes)
        tcPart = 'full';
        leaveOut = '0';
        if strcmp(fitTypes{ff},'V2V3')
            V2V3 = '1'; % '0' = V1-V2, V1-V3; '1' = V1-V2, V1-V3, V2-V3
        else
            V2V3 = '0'; % '0' = V1-V2, V1-V3; '1' = V1-V2, V1-V3, V2-V3
        end
        for ss = 1:length(sessions)
            session_dir = sessions{ss};
            outDir = fullfile(session_dir,'fit_template_scripts',tTypes{tt},fitTypes{ff},volFunc);
            system(['rm -rf ' outDir]);
            mkdir(outDir);
            saveDir = fullfile(session_dir,'pRFs',tTypes{tt},volFunc,'Movie',fitTypes{ff});
            system(['rm -rf ' saveDir]);
            mkdir(saveDir);
            runs = movieRuns{ss};
            create_regress_template_scripts(session_dir,tTypes{tt},outDir,runs,...
                volFunc,saveDir,tcPart,leaveOut,V2V3,logDir);
        end
    end
end
%% Submit the regress scripts (must be run from chead on UPenn cluster) - (coarse)

% e.g.
% cd /data/jet/abock/data/Retinotopy_Templates/AEK/10012014/fit_template_scripts/coarseV2V3size/V2V3/s5.wdrf.tf/
% sh submit_lh_regress.sh
% sh submit_rh_regress.sh
%%

templateType    = 'coarseV2V3size';
fitType         = 'V2V3';
V2V3            = 1;
cluster         = '1'; % run on cluster
tcPart          = 'full';
leaveOut        = '0';

hh = 1:5;
ii = 1:5;
jj = 1:7;
kk = 1:7;
ll = 1:7;

for ss = 1:length(sessions)
    session_dir = sessions{ss};
    runs = eval(movieRuns{ss});
    saveDir = fullfile(session_dir,'pRFs',templateType,volFunc,'Movie',fitType);
    for hhh = 1:length(hemis);
        hemi = hemis{hhh};
        for h = hh
            for i = ii
                for j = jj
                    for k = kk
                        for l = ll
                            if ~exist(fullfile(saveDir,...
                                    [hemi '.' num2str(h) '.' num2str(i) '.' num2str(j) '.' num2str(k) '.' num2str(l) '.varexp.txt']),...
                                    'file');
                                regress_template(session_dir,saveDir,templateType,...
                                    runs,hemi,volFunc,[num2str(h) '.' num2str(i) '.' num2str(j) '.' num2str(k) '.' num2str(l)], ...
                                    cluster,tcPart,leaveOut,V2V3);
                            end
                        end
                    end
                end
            end
        end
    end
end
%% Find the best template (coarse)
hemis = {'lh' 'rh'};
%templateType = 'coarse';
templateType = 'coarseV2V3size';
func = 's5.wdrf.tf';
fitTypes = {'V2V3'};%{'V1' 'V2V3'};
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    disp(session_dir);
    for ff = 1:length(fitTypes)
        fitType = fitTypes{ff};
        disp(fitType);
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            disp(hemi);
            tDir = fullfile(session_dir,'pRFs',templateType,volFunc,'Movie',fitType);
            if strcmp(templateType,'anat')
                vals = load(fullfile(tDir,[hemi '.' templateType '.varexp.txt']));
                if strcmp(fitType,'V1')
                    varexp = nansum(vals(1:2));
                elseif strcmp(fitType,'V2V3')
                    varexp = nansum(vals(1:3));
                end
                disp(varexp);
            else
                [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,fitType);
                disp(params(1));
                disp(sorted_templates(1));
                disp(varexp(1));
            end
        end
    end
end