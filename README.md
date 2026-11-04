# Data Center Heat Analysis

数据中心周边卫星地表温度与粗略回归距离分析。分析对应4590条原始记录。

本仓库包含代码、IDEA运行配置和方法说明。原始CSV、卫星缓存、温度汇总及结果文件保留在本地，不包含在仓库中。下面的目录说明列出完整本地项目所需的数据位置。

## 目录

- `data/US Sites 2026-April (1).csv`：你的原始数据，未修改。
- `work/all_heat.py`：Landsat 数据读取、质量筛选及环带温度提取。包含此前历史窗口处理代码，默认查看和半径重算均不调用它。
- `work/export_all_heat.py`：从已缓存观测汇总温度曲线，并执行严格规则。
- `work/heat_radius_only.py`：整理严格规则结果，按粗略距离档报告。
- `work/relaxed_heat_radius.py`：保留严格结果，并补充宽松首次回归距离。
- `work/impute_all_radii.py`：五折分组验证、中位数插补及来源标注。
- `work/verify_all_heat.py`：核对已有卫星分析的数据和坐标对应关系。
- `work/heat_all_cache/`：已有观测缓存，可离线重新汇总。
- `outputs/heat_all/`：已有卫星观测、环带温度、场景来源及示例热图。
- `outputs/heat_radius_complete/`：完整4590条半径结果，含CSV和已导出的Excel。
- `data_center_analysis.py`：默认显示现有结果概况，不下载数据、不清洗年份、不做DID。

## 在 IDEA 中运行

项目使用 Python。使用项目的 `.venv` 解释器。主要依赖可直接安装：`python -m pip install numpy rasterio pyproj matplotlib`。无需 requirements 文件。
运行配置 `Data Center - Overview` 显示结果概况。
运行配置 `Data Center - Recalculate Radii` 使用已有卫星汇总重新计算严格结果整理、宽松结果和模型插补，无需联网。

在项目终端也可运行：

```powershell
.\.venv\Scripts\python.exe data_center_analysis.py
.\.venv\Scripts\python.exe data_center_analysis.py --recalculate
```

完整重新汇总严格结果：`python work/export_all_heat.py`，然后 `python data_center_analysis.py --recalculate`。
核对已有卫星数据：`python work/verify_all_heat.py`。
重新提取卫星数据的入口为 `python work/all_heat.py`，需要网络访问公共数据服务。参考年份固定为2026；此入口同时保留早期历史窗口逻辑，不是当前默认流程。

## 指标含义

地表温度来自2026年6—8月的Landsat Collection 2 Level-2。各环带取有效陆地温度中位数，背景为6—8公里区域。同一影像先相减得到温差，再汇总日期。距离从CSV坐标点起算，不是园区边界。

严格规则要求12个分析环带各有两次有效观测，近区两个环带均超过背景1°C，从检出位置起后续所有环带均处于背景±1°C，至少连续三个环带。宽松规则允许一次观测，近区可用环带平均温差超过1°C，首次连续两带处于±1°C即接受。

已有1566条卫星半径（522条严格，新增1044条宽松）。其余3024条统一按2.5公里中位数插补，来源分别标记。插补值不证明该站点存在热影响。距离是描述性的地表温度回归背景距离，不是已证实的废热传播边界。

CSV重算会更新结果CSV和方法记录；随附Excel是当前快照，不随重算自动刷新。

来源：https://www.usgs.gov/landsat-missions/landsat-collection-2-surface-temperature
公共数据：https://planetarycomputer.microsoft.com/dataset/landsat-c2-l2


## 半径分布分析

新增脚本 `radius_distribution.py`。在IDEA中直接运行，或在项目终端执行：

```powershell
.\.venv\Scripts\python.exe radius_distribution.py
```

输出在 `outputs/radius_distribution/`：全量与严格规则分布图、按0.5公里档的频数CSV、汇总JSON。
全量4590条含1566条卫星结果和3024条统一2.5公里的中位数插补。图中按来源堆叠，不能把插补峰值解释为实际热影响集中在2.5公里。
严格规则平均值只使用522条识别成功的粗略取整半径，不把未识别记录设为零；当前平均2.4377公里。


## 在新电脑上使用

创建虚拟环境并安装依赖：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install numpy rasterio pyproj matplotlib
```

将本地项目的 `data/`、`outputs/` 和 `work/heat_all_cache/` 复制到相同目录后，才能查看或复现已有结果。代码不包含这些研究数据，也不包含访问令牌。
默认概况脚本读取现有汇总，`--recalculate` 依赖已有卫星温度汇总。它们不直接从原始CSV完成全流程。
严格半径平均值以已识别记录的粗略取整半径计算。所有分布按记录计数，共享坐标或园区的记录没有当作独立设施去重。
