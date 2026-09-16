function status=pmcRunSample(job,options,environment)
% Worker owns one model, its checkpoint and log; large arrays never cross IPC.
started=tic;status=struct('name',job.name,'state','failed','seconds',0,'message','');
logFile=fullfile(options.resultsFolder,'logs',[job.name '.log']);
fid=fopen(logFile,'a');assert(fid>=0,'Cannot open sample log.');
logCleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
oldDir=pwd;dirCleanup=onCleanup(@() cd(oldDir)); %#ok<NASGU>
try
    fprintf(fid,'[%s] START %s\n',datestr(now,30),job.source);fflushLog(fid);
    if options.configureSolver
        restoreEnvironment(environment);
        changeCobraSolver('ibm_cplex','LP',0,-1); % -1 configures globals without returning status.
        global CBT_LP_SOLVER
        assert(strcmp(CBT_LP_SOLVER,'ibm_cplex'),'CPLEX initialization failed.');
    end
    loadStarted=tic;
    if exist(job.cache,'file')
        cached=load(job.cache,'model','sourceKey','info');
        assert(strcmp(cached.sourceKey,job.sourceKey),'Stale lean-model cache.');
        model=cached.model;info=cached.info;clear cached
    else
        [model,info]=pmcReadModel(job.source);
        pmcAtomicSave(job.cache,struct('model',model,'info',info,'sourceKey',job.sourceKey));
    end
    rxns=pmcSelectReactions(model.rxns,options.metList,options.targetRxns);
    fprintf(fid,'LOAD %.3fs, lean %.3f MiB, %d reactions, %d targets, %d coupling rows\n', ...
        toc(loadStarted),info.modelBytes/2^20,info.reactions,numel(rxns),info.couplingRows);fflushLog(fid);
    minFlux=nan(numel(rxns),1);maxFlux=minFlux;completed=false(numel(rxns),1);
    progressFile=[job.output '.progress.mat'];
    if exist(progressFile,'file')
        old=load(progressFile);
        assert(strcmp(old.key,job.key) && isequal(old.rxns,rxns),'Incompatible sample checkpoint.');
        minFlux=old.minFlux;maxFlux=old.maxFlux;completed=old.completed;clear old
    end
    control=struct('PARALLELMODE',1,'THREADS',1,'AUXROOTTHREADS',0);
    if isfinite(options.lpTimeLimit),control.TILIM=options.lpTimeLimit;end
    chunkSize=min(options.reactionChunkSize,max(1,numel(rxns)));
    for first=1:chunkSize:numel(rxns)
        index=first:min(first+chunkSize-1,numel(rxns));
        index=index(~completed(index));
        if isempty(index),continue;end
        fprintf(fid,'[%s] FVA targets %d:%d of %d\n',datestr(now,30),index(1),index(end),numel(rxns));fflushLog(fid);
        t=tic;
        % Capture native CPLEX diagnostics in the sample log, not just Future.Diary.
        solverText=evalc('[lo,hi,objective,ret]=options.fvaFunction(model,options.optPercentage,''max'',''ibm_cplex'',rxns(index),''A'',control,0,[],0);');
        fprintf(fid,'%s',solverText);fflushLog(fid);
        if ret~=0 && options.numericalRetry && strcmp(func2str(options.fvaFunction),'pmcFastFVA') ...
                && contains(solverText,'(status 5)')
            retryControl=control;
            retryControl.NUMERICALEMPHASIS=1;retryControl.SCAIND=1;retryControl.ADVIND=1; % Reuse bases within the newly built LP.
            fprintf(fid,'RETRY: unscaled infeasibility; numerical emphasis and alternate scaling; constraints unchanged.\n');fflushLog(fid);
            solverText=evalc('[lo,hi,objective,ret]=options.fvaFunction(model,options.optPercentage,''max'',''ibm_cplex'',rxns(index),''A'',retryControl,0,[],0);');
            fprintf(fid,'%s',solverText);fflushLog(fid);
        end
        cd(oldDir);
        assert(isscalar(ret) && ret==0,'PMC:FVAStatus', ...
            'FVA return code %g for target positions %d:%d. See native solver diagnostics in sample log.',ret,index(1),index(end));
        assert(numel(lo)==numel(index) && numel(hi)==numel(index),'PMC:FVASize','Wrong FVA output size.');
        assert(all(isfinite(lo)) && all(isfinite(hi)) && all(lo(:)<=hi(:)+1e-6), ...
            'PMC:FVAValues','Invalid/unfinished FVA values; not saved as complete.');
        minFlux(index)=lo;maxFlux(index)=hi;completed(index)=true;
        fprintf(fid,'DONE chunk %.3fs\n',toc(t));fflushLog(fid);
        if isfinite(options.reactionChunkSize)
            pmcAtomicSave(progressFile,struct('key',job.key,'rxns',{rxns}, ...
                'minFlux',minFlux,'maxFlux',maxFlux,'completed',completed));
        end
    end
    result=struct('key',job.key,'sourceKey',job.sourceKey,'modelFile',job.name, ...
        'source',job.source,'rxns',{rxns},'minFlux',minFlux,'maxFlux',maxFlux, ...
        'state','complete','optPercentage',options.optPercentage,'info',info, ...
        'seconds',toc(started),'origin','computed','solverReturnCode',0,'schemaVersion',1);
    pmcAtomicSave(job.output,struct('result',result));
    if exist(progressFile,'file'),delete(progressFile);end
    status.state='complete';fprintf(fid,'[%s] SAVED %s\n',datestr(now,30),job.output);
catch err
    status.seconds=toc(started);
    status.message=getReport(err,'extended','hyperlinks','off');
    fprintf(fid,'[%s] FAILED\n%s\n',datestr(now,30),status.message);
    pmcAtomicSave([job.output '.failure.mat'],struct('status',status,'key',job.key));
end
status.seconds=toc(started);
end
function fflushLog(fid)
% fseek on a regular output file flushes MATLAB's buffered stream in R2019b.
fseek(fid,0,'cof');
end
