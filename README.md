# Group-meeting
## Introduction

组会代码整理，将每次组会报告中编写的程序进行整理。

### 2023.3.16
Unbounded quasiperiodic models and anomalous mobility edges.\
Perturbation calculation for Multifrctal states.\
The origin of MBC region and Ergodic breaking. (Python code and need to be future organized.)

### 2023.6.11
strong disorder renormalization group\
ground state in Heisenberg model\
RSRG-X calculate the excited states\
Avalanche in MBL can also use the SDRG.\ 

## 实际代码入口与复现范围

| 路径 | 当前入口/作用 |
| --- | --- |
| `20230316/Unbounded/MarylandModel.m` | Maryland 模型的谱与本征态 |
| `20230316/Unbounded/UnboundedQuasiperiodicModels.m` | 参数扫描与能量/IPR 图 |
| `20230316/MultifractalPerturbation/ChangeLambda.m` | 默认运行 zeropoints 分析 |
| `20230316/MultifractalPerturbation/PerturbationMain.m` | 默认运行 PerturbationAnalytic4 分析 |
| `20230316/组会20230316.pdf` | 对应组会报告 |
| `Strong Disorder RG/20230608RSRG.pdf` | SDRG 报告；此目录当前没有完整算法代码 |

使用支持脚本内局部函数的 MATLAB（R2016b 或更新版本）。
从仓库根目录运行时，可先添加对应实验目录：

```matlab
experiment = fullfile(pwd,'20230316','Unbounded');
addpath(experiment);
rng(0);  % 可选，固定随机相位
run(fullfile(experiment,'MarylandModel.m'));
```

两个多分形脚本共享其目录中的 `Ipr.m`；无界模型脚本复用已存在的
同目录 `Ipr.m`。计算公式保持原样。跨仓库的历史函数副本仍然保留，
运行此仓库不需要另一个 checkout。部分脚本默认矩阵较大，
应在运行前检查 L 和参数扫描范围。标题中的 MBC/Python 待整理说明
是历史研究计划，不代表当前目录已包含这些程序。`test` 占位文件
也不是测试套件。

```matlab
results = runtests('tests');
assert(all([results.Passed]));
```

测试只覆盖 IPR 函数与本地路径解析；不验证组会报告的物理结论。

