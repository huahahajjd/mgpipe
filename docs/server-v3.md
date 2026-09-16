# 标的反应 + 连续补位（v3，2026-09-15）

已对全部 418 个原始模型读取矩阵维度：352 个 MAT v5/v7 文件通过流式解压读取元数据，66 个 v7.3 文件通过 MATLAB HDF5 接口读取稀疏矩阵元数据。反应数按 S 的列数、耦合约束数按 C 的行数统计，不按文件大小排序。

唯一同时达到两个最小值的标的为 `ERR3472679`（原排序第 12 个）：59,405 个反应变量、83,915 条耦合约束、53,136 条代谢平衡行。再次实际加载模型核验维度后，得到并保存了 7,928 个完整 IEX 反应 ID。

## 分析范围

每个样本只对“本样本 IEX 反应与标的 IEX 反应 ID 的交集”进行 FVA。使用完整 ID（包含菌株前缀），不是仅按代谢物名称匹配。不在标的中的菌株反应不再计算。本样本不存在的标的反应不会凭空添加；选中目标为空时保存空结果。

所有样本仍保留自己的完整 S、C、d、dsense、边界和目标。此处减少的是要优化的目标反应数，不是把所有样本改造成标的的网络。科学结果的覆盖范围缩小了，不能将新表解释为该样本全部菌株的贡献。

## 调度

正式入口是 `run_predict_microbe_contributions_v3`。固定 8 个进程、滚动连续补位：任务结束就领取后续任务，不等待 8 个或 100 个样本组成的批次。入口没有 batchSize 参数。兼容函数仍接受旧批次选项，但 v3 固定使用 rolling。

每个样本完成/失败后立即更新 run_summary.mat；失败不阻止其他样本继续，队列会处理到最后。仍采用按反应分块（每块至多 1000 个目标）保存断点，这不是样本批次。

CPLEX 状态 5 会记录原始消息，并对该块使用数值精度强调、替代缩放重新构建 LP 求解一次（保留块内基复用）。网络约束及容差不放宽；重试仍不满足最优状态就保留失败及断点，不填零、不标为完成。其他类型错误不会无限重试。成功结果全部完成后才导出全体 CSV；有失败则保留逐样本结果和失败清单，重跑时只调度未完成样本。

## 数据保存与旧结果迁移

旧任务在 2026-09-15 13:54:27（服务器时间）停止，退出码 -2。停止前备份 135 个 MAT 文件及样本日志至：

`/mnt/nfs/wangchao/my_course/untreated_MS/replace_data/step4_modeling/predict_microbe_contributions/stop_checkpoints/20260915_135401`

新结果目录为 `run_reference_rolling_20260915`。旧目录 `run_refactored_8` 保留，用于验证同一源模型、原目标比例和旧 checkpoint key 后迁移：完整结果取目标交集；部分结果仅迁移完成标志为 true 的目标。新范围已全部算完的样本直接标为完成，其余保留部分进度。原恢复结果继续保留 origin=recovered、未知原求解状态码的说明。

缓存仍共用 `lean_model_cache_v1`。原始样本模型不覆盖。

## 运行

```matlab
base='/mnt/nfs/wangchao/my_course/untreated_MS/replace_data/step4_modeling/predict_microbe_contributions';
addpath(base);
summary=run_predict_microbe_contributions_v3();
```

同一命令用于续算。不要同时启动两个协调器写同一目录；后台启动记录由 `latest_production_launch.json` 定位，退出状态由对应启动目录的 `exit_status.json` 记录。

完整清单在 `reference_targets_20260915.mat`，规模统计在开发目录的 `model_dimensions.json`。工作状态以最新工作断点文档和正式启动日志为准。
