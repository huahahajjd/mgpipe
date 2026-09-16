# mgpipe — 微生物贡献预测 v3

基于 COBRA Toolbox 的 `predictMicrobeContributions` 重构，针对多个个性化群落模型进行可续算的 FVA。核心源码与 2026-09-15 部署的 v3 一致。

## 特性与分析范围

- 默认 8 个 MATLAB worker，rolling 连续补位，复用并行池（IdleTimeout=Inf），无样本批次等待。
- 使用标的模型的完整 IEX 反应 ID（包含菌株前缀），每个样本仅优化与标的的交集。保留完整代谢网络、边界及耦合约束；结果仅覆盖选中的反应。
- 精简模型缓存，减少注释等非求解字段的读取和内存开销；不覆盖原始模型。
- 样本独立保存，支持反应块断点、已有结果迁移、内存准入及失败隔离。
- CPLEX 单线程、基复用和每个 LP 的最优状态检查。状态 5 时用数值强调与替代缩放重试一次，失败仍如实记录，不放宽约束。

本项目 418 个模型中，共同最小标的为 ERR3472679（59,405 个反应、83,915 条耦合约束），对应 7,928 个 IEX ID。其他数据集应重新选择标的。减少目标数不等于缩小单次 LP 的规模。

## 环境与安装

已验证环境：Linux x86-64、MATLAB R2019b、Parallel Computing Toolbox、COBRA Toolbox、IBM CPLEX 12.9，以及 MATLAB MEX 支持的 C 编译器。CPLEX 和 MATLAB 需自行安装；本仓库不包含其二进制或许可证。

```matlab
repo = '/path/to/mgpipe';
cobraRoot = '/path/to/cobratoolbox';
cplexRoot = '/path/to/CPLEX_Studio';
addpath(fullfile(repo,'src'));
pmcConfigureCobra(cobraRoot,cplexRoot);
addpath(fullfile(repo,'src'),'-begin'); % 优先使用本仓库版本
build_pmc_backend(cplexRoot);           % Linux/CPLEX 12.9 后端
which predictMicrobeContributions
```

## 运行与续算

以下示例使用已选择的标的模型生成目标列表。标的应同时具有最少反应变量和耦合约束；若不存在共同最小值，需要明确选择标准。

```matlab
modelDir = '/path/to/models';
[referenceModel,~] = pmcReadModel('/path/to/reference_model.mat');
targetRxns = pmcSelectReactions(referenceModel.rxns,{});
assert(~isempty(targetRxns),'Reference target list is empty');
[~,~,~,summary] = predictMicrobeContributions(modelDir, ...
    'resultsFolder','/path/to/results', ...
    'cacheFolder','/path/to/cache', ...
    'targetRxns',targetRxns, ...
    'scheduling','rolling','numWorkers',8, ...
    'reactionChunkSize',1000,'numericalRetry',true, ...
    'failOnError',false,'returnTables',false);
```

重复同一命令即可续算。先加 `'dryRun',true` 可预检；预检会创建目录、清单并执行启用的迁移，但不求解。`targetRxns={}` 表示全部 IEX，不表示空集合。省略 `reactionChunkSize` 时函数默认 Inf；以上入口显式设为 1000。`failOnError` 函数默认 true，以上显式设为 false 以处理完整队列。

默认每个活跃样本预算 8 GiB，额外保留 8 GiB；可用内存不足时活跃任务数会降至 8 以下。池的无限空闲时间不防止主进程退出。不要让两个协调器写同一结果目录，也不要在计算中修改源模型。异常退出留下锁时，先确认旧进程和 worker 均已停止。

结果在 `samples/`，样本日志在 `logs/`，总体状态在 `run_summary.mat`。反应块断点为 `.progress.mat`；`.failure.mat` 为历史诊断，最终状态以正式结果为准。全部样本成功后才导出总体 CSV。失败样本会在下次续算时重试。

## 目录与复现

- `src/`：核心函数及可自行编译的 CPLEX MEX 源码。
- `tests/`：数值、调度、I/O、断点和迁移测试；部分集成测试含原服务器路径，需修改后运行。
- `diagnostics/`：原项目 MAT v5/v7、v7.3 模型规模统计脚本，含项目路径和数量假设。
- `examples/server/`：原服务器 v3 入口、参考目标准备、安装与后台启动脚本，原样归档。**不能直接当作通用安装器执行**：安装器依赖旧版哈希清单，启动器依赖部署清单，准备脚本要求 418 个模型。移植时必须修改路径和前提条件。
- `docs/server-v3.md`：原部署说明；其中运行目录、备份及数据文件仅在原服务器存在。
- `docs/source-manifest.json`：核心源码 SHA-256，用于对照已部署版本。

## 验证与已知限制

服务器部署前已通过真实 CPLEX 小模型数值、完整/部分结果迁移、rolling 失败补位、I/O 和断点测试。本次仓库整理不代表重新执行 MATLAB 测试。初始化环境并添加 `tests/` 后可运行 `test_pmc_io`、`test_pmc_scheduler`、`test_pmc_reference`、`test_pmc_resume`；测试会创建并行池并需要相应资源。

生产运行仍有少数样本在数值重试后返回 CPLEX 状态 5，v3 会保留失败记录并继续队列，未宣称解决所有数值问题。原内存恢复结果的求解状态未知，迁移时保留来源标记。完整模型、结果、恢复文件、缓存及运行日志未收录。

## 来源与许可

`src/pmcCplexFVA.c` 改编自 fastFVA/CPLEXINT，保留作者、版权及 LGPL-2.1-or-later 声明，修改说明见文件头；许可文本见 `LICENSES/`。COBRA Toolbox、MATLAB、CPLEX 为外部依赖。本仓库未另行为其余文件授予统一开源许可证。
