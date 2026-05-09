# hn_work_flow 使用说明

## 1. 项目简介

`hn_work_flow` 是一个面向分子动力学任务的脚本化工作流，核心目标是把下面三件事串起来：

1. 自动建模并生成一批可用于后续统计的初始结构。
2. 批量创建 `case1 ... caseN`，自动生成作业脚本并提交扫描电压/充电类 MD 任务。
3. 对多个 case 的结果做统一后处理与平均。

当前仓库的主入口是 `md_flow.sh`，实际可用的主模块为：

- `frames`：建模并提取初始结构
- `run`：批量建 case 并提交任务
- `analyze`：对 `case1 ... caseN` 做后处理
- `d_data`：历史备份模块，依赖 Matlab，通常可以忽略

这个工作流主要围绕 GROMACS + shell 脚本组织，局部使用 C++ 程序做文件检查、建模辅助和后处理。

## 2. 仓库结构

```text
.
├─ md_flow.sh              # 总入口，读取 INPUT 后分发到各模块
├─ INPUT                   # 唯一配置入口
├─ INPUT_guide             # INPUT 参数说明
├─ build_model/            # 建模与初始结构提取
├─ run_md/                 # 建 case、生成作业脚本、提交任务、rerun
├─ analyze/                # 后处理分析框架
├─ deal_data/              # 旧数据处理模块，主要作为备份
├─ example/                # 各模块示例
└─ bin/                    # 编译工具/可执行文件
```

`example/` 中已经给了完整示例：

- `example/channel/set_density`：电极体系，直接指定 `set_density`
- `example/channel/get_density_from_bulk`：电极体系，通过 bulk NPT 获取目标密度
- `example/bulk`：体相建模
- `example/slit`：狭缝孔建模
- `example/run_md`：批量运行任务
- `example/analysis`：分析模块

## 3. 依赖与运行环境

工作流本身无需安装，但运行环境至少应具备：

- `bash`
- `gromacs`
- `g++`
- `packmol`
- 对应集群的提交命令
  - `eninstein`：`qsub`
  - `4090`：`sbatch`
  - `wuchao`：`dsub`
  - `jiaocha`：`sbatch`
- `Matlab`，仅 `d_data` 模块需要

首次使用建议：

```bash
chmod +x md_flow.sh
```

之后可把仓库路径加入 `PATH`，即可在任意位置调用

## 4. 基本使用方式

```bash
./md_flow.sh -h
./md_flow.sh [module_name] [additional_arguments]
```

### 4.1 可用命令

```bash
./md_flow.sh frames
./md_flow.sh frames bulk
./md_flow.sh frames slit
./md_flow.sh frames clear
./md_flow.sh frames clear --dry-run

./md_flow.sh run
./md_flow.sh run rerun

./md_flow.sh analyze eleQ
./md_flow.sh analyze onfly
./md_flow.sh analyze mdheat
./md_flow.sh analyze MD_PP
./md_flow.sh analyze MD_PP rerun
```

说明：

- `frames` 默认是电极体系 NVT 建模。
- `frames bulk` 是体相建模。
- `frames slit` 是狭缝孔建模。
- `frames clear` 会清理 `frames/build_model` 生成的目录、断点文件、GROMACS 输出和临时文件；可先用 `frames clear --dry-run` 预览。
- `run rerun` 会走 `run_md/rerun.sh`。
- `analyze` 当前脚本里真正支持的分支名是 `eleQ`、`onfly`、`mdheat`、`MD_PP`。
- 如果你需要新增了自己的分析分支，可以在 `analyze/analyze.sh` 中添加。

## 5. 整体工作流

### 5.1 建模模块：`frames`

适用范围：

- 电极体系 NVT 建模
- 体相 NPT 建模
- 狭缝孔 NVT 建模（包括新建电极pdb文件）

默认流程如下：

1. 根据 `set_density` 调整体系密度，循环修改指定分子数量并做 NVT 模拟检查密度，模拟步数对应 `nsteps_den`。
2. 若未设置 `set_density`，工作流会先进行 bulk NPT（nsteps_bulk），自动得到平衡密度，再把该密度作为目标密度继续后续流程。
3. 调密度完成后继续做一段平衡模拟，对应 `nsteps_equ`。
4. 最后再跑一小段轨迹，对应 `nsteps_ini`，从中筛选出最接近设定温度的结构。
5. 输出 `num_intialframes` 个初始结构，命名为 `frame1.gro ... frameN.gro`。

实际输出目录：

- 电极体系/狭缝孔：`./{Temperature}k`
- 体相：`./{Temperature}k_bulk`

补充说明：

- 体相建模时，程序会自动做 bulk NPT，然后跳过电极体系的调密度步骤；此时 `set_density` 无效。
- 狭缝孔建模会比普通电极建模多出新建电极步骤。
- `build_model/run.sh` 和 `build_model/build_slit.sh` 已支持断点续跑，分别通过 `.step.cpt` 和 `.build_slit.step.cpt` 记录进度。

### 5.2 任务提交模块：`run`

适用范围：

- 适配使用 CPM 充电的扫描电压任务
- 其他类型任务也能用，但可能会生成一些冗余文件或文件夹
- 如需改写生成扫描电压文件的算法（目前为斜坡电压），可修改 `run_md/source/generate_scan.cpp`

模块行为：

1. 复制初始结构、矩阵文件和基础文件到目标目录。
2. 自动生成 `case1 ... caseN`。
3. 生成并提交集群作业脚本。

默认输出目录：

```text
./charging/{Temperature}k/{V}V/{ic}ps/case1...caseN
```

`rerun` 模式下，会在已有 case 中生成 `rerun_case` 并重新投递。

### 5.3 分析模块：`analyze`

适用范围：

- 电极电荷统计
- onfly 二进制输出处理
  - 提供三种版本的 onfly 处理程序
- 体系产热计算
- 用户自定义后处理，内置范围包括：
  - 三维数密度动态分析
  - 分子偶极矩与z轴夹角的分析
  - 动态电场分析
  - 沿 z 方向统计分子电荷、偶极、四极矩，多极展开分析

对应于：

- `eleQ`：平均指定 case 范围的 `CPM_electrodeCharge.dat`
- `onfly`：平均指定范围的 onfly 输出
- `mdheat`：体系产热相关计算
- `MD_PP`：调用自定义 C++ 后处理程序

这个模块本质上是一个开放框架。推荐的组织方式是：

- 在 `analyze/analyze.sh` 中只做分支分发
- 下一层 shell 脚本负责参数整理
- 真正的计算逻辑放在 `analyze/source/module_C++` 里的 C++ 程序中

必须遵守的目录约定：

- 模拟目录需要满足 `case1`、`case2`、...、`caseN`，且 mdrun 时使用了 -deffnm 参数。如果使用本仓库的 `run` 模块生成，这个约定会自动满足。
- `INPUT` 需要放在这些 `case*` 文件夹的同级目录

`MD_PP` 的工作方式：

- 通过 `INPUT` 中的 `analyze_cpp` 指定要运行的 C++ 文件名
- 对应源码需放在 `./analyze/source/module_C++/`
- 会自动编译到 `./bin/`
- 内置了 `-f -s -n -o` 等输入输出参数
- 轨迹和 `tpr` 文件名由 `DEFFNM_analyze` 与 `xtcORtrr` 控制

等价命令形式大致为：

```bash
echo ${analyze_mol} | <analyze_cpp> \
  -f ./case1/${DEFFNM_analyze}.${xtcORtrr} \
  -s ./case1/${DEFFNM_analyze}.tpr \
  -n ./case1/index.ndx \
  -o ./deal_data/${analyze_cpp}/1${analyze_cpp}.xvg \
  ${analysis_extra_command}
```

如果是 `rerun` 分析，则分析目录会切换到 `./case1/rerun_case/`。

仓库当前已提供的 `MD_PP` 示例程序包括：

- `3DdensityZl_dynamic2020`
- `angleOF_dipole_Z`
- `dynamic_force_bin2020`

## 6. 各模块所需文件

### 6.1 `frames` 所需文件

1. `INPUT`
2. 拓扑文件
   - `TOP.top`
   - `bulk_top.top`
   - 其中要包含所需 `itp`
3. `mdp` 文件
   - `MDP.mdp`
   - `mini_MDP.mdp`
4. `packmol` 输入文件
   - `packmol.inp`
   - `bulk_pa.inp`
5. 所有需要的 `pdb`

重要注意事项：

- `INPUT` 里的文件名都写“前缀”，不要带扩展名。
  - 正确示范：`TOP = topol`
  - 错误示范：`TOP = topol.top`
- `packmol` 中除电极外的 `pdb` 文件名必须与 `top` 中分子名保持一致。已通过 `build_model/check_file.sh` 来检查这个一致性。
- 如果不设置 `set_density`，需要额外准备 bulk 建模所需的 `top` 和 `inp` 文件。
- 如果只做 bulk 建模，`MDP` 中仍建议保留一行带注释的 `ewald-geometry = 3dc`，删掉或许会出现未知错误。
- 如果不指定 `set_density`，请在 `MDP` 中写好压力耦合相关设置，并在模拟开始前保持关闭状态。

狭缝孔额外说明：

- 狭缝孔建模不需要电极 `pdb`。
- 需要在 `INPUT` 里补充狭缝孔参数。
- 电极分子名固定使用 `EL`、`ER`、`GRA`。
- 单侧电极只包含 `EL` 和 `GRA`，所以初始 `mdp` 的冻结组只需要覆盖这些组。初始mdp不要添加ER！否则会报错！
- 首次使用狭缝孔建模时强烈建议打开 `is_auto_top=1` ，会自动生成电极相关 `top/inp`。
- 如果不是首次使用，务必关闭 `is_auto_top`，避免覆盖你已经整理好的输入文件。

### 6.2 `run` 所需文件

1. `INPUT`
2. `initialframes/`
   - 放 `frame*.gro`
3. `matrixdata/`
   - 放 CPM 相关矩阵文件和控制文件
   - 如果不做 CPM，建议仍提供一个空文件夹，避免未知错误
4. `basicfile/`
   - 放力场、`index.ndx`、`grompp.mdp` 等基础文件
   - 放拓扑文件

参考示例：`./example/run_md`

### 6.3 `analyze` 所需文件

1. `INPUT`
2. `case1 ... caseN`
3. 每个 `case` 下至少包含：
   - `${DEFFNM_analyze}.tpr`
   - `${DEFFNM_analyze}.xtc` 或 `${DEFFNM_analyze}.trr`
   - `index.ndx`

参考示例：`./example/analysis`

### 6.4 `d_data` 模块

`deal_data` 当前主要作为历史备份使用，其中 JE 相关脚本还有特定机器依赖。大多数使用场景下可以忽略这一部分。

## 7. INPUT 参数说明

这部分按 `INPUT_guide` 的风格，对最常用参数做汇总说明。

### 7.1 `frames` 部分

#### `# Simulation box dimensions`

- `xbox ybox zbox`
  - 普通电极体系的盒子尺寸，对应最终 `.gro` 最后一行
  - 因为 `packmol` 生成的 `pdb` 不带盒子信息，所以这里必须显式指定
- `bulk_xbox bulk_ybox bulk_zbox`
  - bulk 体系盒子尺寸
  - 即使是狭缝孔建模，只要你不设置 `set_density`，仍然需要这组 bulk 尺寸参数
- 狭缝孔建模时，电极文件会自动带正确盒子尺寸，此时 `xbox ybox zbox` 可忽略

#### `# input file prefixes (without extensions)`

- `TOP`
- `bulk_top`
- `MDP`
- `mini_MDP`
- `packmol`
- `bulk_pa`

这些都只写前缀，不写后缀名。

#### `# gromacs settings`

- `GMXRC`
  - 建模模块使用的 GROMACS 环境脚本路径
- `NPOS`
  - 建模时 OpenMP 线程数
- `GPU`
  - 是否启用 GPU
- `maxWarn`
  - `gmx grompp` 的 `-maxwarn`
- `nsteps_bulk`
  - bulk NPT 获取密度的步数
- `nsteps_den`
  - 电极体系调密度时单轮 NVT 步数
- `nsteps_equ`
  - 平衡阶段步数
- `nsteps_ini`
  - 提取初始帧前的最后一段运行步数
- `XOUT_FRAMES`
  - 输出轨迹的时间间隔，太小会导致 I/O 负担很重
- `num_intialframes`
  - 最终导出的初始结构数量

#### `# Density analysis parameters`

- `NZ`
  - `gmx density` 中 `-nbin`
- `BE_TIME`
  - 密度分析起始时间
- `BE_Z` / `END_Z`
  - 分析密度时采用的 z 方向范围
- `error`
  - 允许的密度误差百分比
- `max_iter`
  - 调密度最大迭代次数
- `onlyUP`
  - 是否只接受高于目标值一侧的设定

#### `# slit parameters`

- `LLX`
  - 电极厚度
- `LLZ`
  - 狭缝长度
- `dslit`
  - 狭缝宽度
- `len_bulk`
  - 体相区域长度，通常建议 8-10 左右
- `is_auto_top`
  - 是否自动生成狭缝孔的拓扑与 `inp`
- `MOL`
  - 需要插入的分子名
- `MOL_num`
  - 狭缝孔体系中各分子数量
- `MOL_num_bulk`
  - 不设置 `set_density` 时，bulk NPT 阶段使用的分子数量

#### `# set_density`

- `set_density`
  - 电极体系 NVT 调密度的目标值，单位 `kg/m^3`
  - 可以不写
  - 如果不写，程序会先跑 bulk NPT 自动求得平衡密度

#### `# modify molecule numbers`

- `MOL_name`
  - 调密度时允许修改数量的分子名
- `FIXED_MOL`
  - 需要保持固定比例关系的分子
- `IS_FORCE_SCALE`
  - 在缩放后分子数变化过小的情况下，是否强制至少改动 1 个分子

这里有两个非常重要的经验性约束：

1. `FIXED_MOL` 是敏感参数，不要随意开启。
   - 原因是比例关系是离散整数约束，不是连续变量。
   - 如果初始配比不好，可能永远无法收敛到目标密度。
2. `IS_FORCE_SCALE` 也不要无脑开启。
   - 尤其和 `FIXED_MOL` 同时使用时，很容易破坏本来的比例约束并导致无法收敛。

### 7.2 `run` 部分

#### `# input file folder prefixes`

- `GMX_run`
  - 实际执行 `mdrun` 的 gromacs 可执行文件
- `initialframes`
  - 初始结构目录
- `matrixdata`
  - CPM 矩阵目录
- `basicfile`
  - 力场、拓扑、`index.ndx`、`mdp` 等基础文件目录

#### `# fold settings`
 注意这里的参数不会影响真实的模拟设置，只会反映在生成的目录结构和文件命名中，或许未来会和扫描电压参数合并。
- `V`
  - 电压
- `ic`
  - 扫描时间或上升时间参数，目录名里会用到
- `Temperature`
  - 温度，目录名里会用到

#### `# job settings`

- `START` / `END`
  - case 编号范围
- `num`
  - 每个作业打包提交多少个 case
- `queue`
  - 集群队列名
- `mode`
  - `default` 或 `append`，`append`会开启 md 的续跑模式
- `DT`
  - 时间步长
- `DEFFNM`
  - 模拟输出文件前缀
- `server_machine`
  - 当前支持 `eninstein`、`4090`、`wuchao`、`jiaocha`
- `server_core`
  - 节点核数设置

说明：

- 不同集群模板分别在 `run_md/*.job`。
- `sub_job.sh` 会按 `server_machine` 选择模板并写入参数。
- 最终 `-ntomp`、节点资源等配置仍然与 job 模板联动，改资源配置时建议同时检查模板文件。

#### `# job control settings`

- `ismakecase`
  - 是否生成 case 目录
- `ismakejob`
  - 是否生成并提交作业

#### `# scan settings`

- `TIME`
  - 扫描电压文件总时长
- `SKIPTIME_PS`
  - 前段不加电的时长
- `TAO`
  - 电压上升时间，可写成一个或多个值

#### `# onfly settings`

- `MOL_name_run`
  - onfly 中选择的分子
- `LOW_ONFLY`
  - 下边界
- `UP_ONFLY`
  - 上边界
- `NBIN_ONFLY`
  - 分箱数
- `MODE_ONFLY`
  - onfly 模式
- `ONFLY_FLAGS`
  - onfly 输出内容
  - `0(x)` `1(v)` `2(f)` `3(xv)` `4(xf)` `5(vf)` `6(xvf)`

#### `# rerun settings`

- `isRerun`
  - 是否启用 rerun 逻辑
- `rerun_start` / `rerun_end`
  - rerun 的 case 范围

### 7.3 `analysis` 部分

#### `# onfly analysis settings`

- `analysis_begin_t`
  - onfly 分析起始时间
- `analysis_end_t`
  - onfly 分析结束时间，`0` 表示直到末尾
- `analysis_ONFLY_isD`
  - 是否按 dynamic 模式处理
- `analysis_ONFLY_in`
  - onfly 输入文件前缀
- `MODE_ONFLY_analysis`
  - onfly 分析模式

说明：

- 文档语义上这里应由 `MODE_ONFLY_analysis` 控制。
- 但当前 `analyze/source/get_onfly.sh` 实际使用的是 `MODE_ONFLY`。如果你的分析结果异常，先检查这两个参数是否一致。

#### `# analysis case settings`

- `analyze_begin_case` / `analyze_end_case`
  - 分析的 case 范围
- `analyze_T`
  - 输出命名中使用的温度
- `analyze_V`
  - 输出命名中使用的电压
- `analyze_tau`
  - 输出命名中使用的时间参数
- `analyze_mol`
  - 要分析的分子名
- `DEFFNM_analyze`
  - 轨迹与 `tpr` 的公共前缀

补充：

- 如果 `index.ndx` 中有重名分子，`MD_PP` 分支会优先取第一个。
- 如果因为重名导致失败，脚本会尝试执行 `fix_duplicate_molname.sh` 后再次运行。

#### `# Analysis program settings (MD_PP)`

- `analyze_core`
  - 并发执行的分析任务数
- `analyze_cpp`
  - 自定义后处理程序名
- `build_gmx`
  - 编译该程序所需的 GROMACS 环境脚本
- `xtcORtrr`
  - 轨迹格式，填 `xtc` 或 `trr`
- `analysis_extra_command`
  - 额外传给后处理程序的参数
- `analyze_is_skip_first_line`
  - 少数程序平均时是否跳过首行，通常保持 `0`

## 8. 日志与输出

主日志默认都放在当前工作目录的 `result/` 下：

- `result/b_model.log`：建模模块日志
- `result/run_md.log`：任务生成/提交日志
- `result/analyze.log`：分析模块日志
- `setting.log`：`INPUT` 读取后的导出记录

常见输出目录：

- `./{Temperature}k` / `./{Temperature}k_bulk`：初始结构
- `./charging/.../caseN`：运行 case
- `./deal_data/`：分析输出

## 9. rerun 的使用

`run` 和 `analyze` 都支持 rerun。

前提：

- 当前目录下已经存在 `case1 ... caseN`

常用命令：

```bash
./md_flow.sh run rerun
./md_flow.sh analyze MD_PP rerun
```

这时：

- `run` 会生成或提交 `rerun_case`
- `analyze` 会改为读取 `./caseN/rerun_case/` 下的文件

## 10. 常见使用建议

1. 先从 `example/` 中最接近你体系的示例目录开始改。
2. 先确认 `INPUT` 里的文件名前缀写法正确，再跑 `frames`。
3. 电极体系调密度时，非必要少用 `FIXED_MOL` 和 `IS_FORCE_SCALE`。
4. 狭缝孔第一次建模时开 `is_auto_top`，后续复用时关掉。
5. 分析模块务必保证目录命名严格是 `case1 ... caseN`。
6. 自定义后处理时，优先遵循当前框架：`analyze.sh` 分发，下一层 shell 整理参数，C++ 程序做计算。

## 11. 一句话上手顺序

如果你是第一次用这个仓库，可以使用以下顺序：

1. 复制一个 `example/` 目录到新的工作目录。
2. 修改 `INPUT` 和对应的 `top/mdp/inp/pdb` 文件。
3. 运行 `./md_flow.sh frames` 或 `./md_flow.sh frames bulk/slit`。
4. 确认得到 `frame*.gro` 后，再准备 `initialframes/ matrixdata/ basicfile/` 运行 `./md_flow.sh run`。
5. 任务结束后，在 `case1 ... caseN` 同级目录放好 `INPUT`，再运行 `./md_flow.sh analyze ...`。

## 12. `INPUT_guide` 补充说明

下面这些提示来自 `INPUT_guide`，因为更偏参数级说明，所以单独补在末尾，避免打乱前面的主体结构。

### 12.1 `frames / Density analysis parameters`

- `NZ` 和 `BE_TIME` 在实际使用上可直接对应到 `gmx density -b {BE_TIME} -nbin {NZ}`。
- `error`、`max_iter`、`BE_Z`、`END_Z` 一起决定调密度时的误差阈值、最大迭代次数和密度统计区间。

### 12.2 `run / job settings`

- 运行任务时，`mdrun` 的命令形态可理解为：

```bash
gmx mdrun -s {DEFFNM}.tpr -deffnm {DEFFNM} []/[-cpi {DEFFNM}.cpt -append] -v -ntmpi 1 -ntomp {硬编码}
```

- 其中 `-ntomp` 的具体写入位置见 `./run_md/sub_job.sh` 及对应的 job 模板。

### 12.3 `run / onfly settings`

- `ONFLY_DENSITY3D` 的环境变量格式示例为：

```bash
export ONFLY_DENSITY3D="-low [] -up [] -nbin [] -n index.ndx -sel [Na CLO aACN]"
```

- `MODE_ONFLY` 一般只有在使用改写版 onfly 时才需要额外关注，默认建议保持为 `3`。

### 12.4 `analysis / onfly analysis settings`

- 分析阶段涉及 onfly 时，`MODE_ONFLY` 需要和你实际运行任务时使用的 onfly 模式保持一致。
- `analyze_mol` 若对应多个index中重名分子，脚本默认优先取第一个分子名。

### 12.5 保留 `INPUT_guide`

- 当前 README 已经合并了 `INPUT_guide` 的主体信息，但为了保守起见，本轮仍保留 `INPUT_guide` 文件，不做删除。
- 如果后续确认 README 与 `INPUT_guide` 完全等价，再删除会更稳妥。
