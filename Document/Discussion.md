还是我们最开始的讨论的论文核心工作，方法学创新，写作技巧、方式要参考E:\CausalSDMs\references\CISO__Species_distribution_modelling_Conditioned_on_Incomplete_Species_Observations.md 这篇论文，方法学要突出我的核心工作是什么？我做了别人没有做！

CAST最后的目标是一个通用的物种分布建模（（像CISO那样通过什么东西解决了什么问题形成了建模方法），不单单针对鱼类，鲫鱼实验仅仅是测试用例；

CAST 的创新远不止变量筛选——它还在建模端做了创新（DAG-Net / CPF），这部分才是最有论文区分度的。

CAST 应该是：
一个完整的因果驱动SDM框架，串接为：
数据 → 因果结构学习(DAG) → 因果效应估计(ATE/CATE) → 因果变量筛选 → 因果结构化建模(DAG-Net)


SDM变量选择中的因果盲区。 传统SDM将所有环境变量视为"平等的输入"，通过统计关联（VIF共线性、LASSO正则化）筛选变量。但统计关联 ≠ 因果关系——被保留的变量可能是混杂代理，被剔除的可能是真正驱动因子。当模型外推到新环境时，这些虚假关联就会失效


方案1：因果先验森林（Causal Prior Forest, CPF）⭐⭐⭐⭐
核心思想：标准Random Forest在每棵树的每次分裂时等概率随机选变量。CPF改为用因果权重加权选择——DAG出度高、ATE显著的变量被优先用于分裂。

标准RF:  P(选Xⱼ做分裂) = 1/p  (均匀)
CPF:     P(选Xⱼ做分裂) ∝ α·OutDeg(Xⱼ) + β·|ATE(Xⱼ)| + γ  
         因果驱动因子更容易被选为分裂变量
为什么好：

ranger R包原生支持 split.select.weights 参数，一行代码就能实现
理论依据清晰：因果驱动变量应该比乘客变量获得更多建模资源
与传统RF完美对比：同样的数据、同样的算法框架、唯一区别是分裂权重来自因果推断
难度：★☆☆☆☆（几乎不需要写新代码）


方案4：DAG结构化神经网络（DAG-Net）⭐⭐⭐⭐⭐
核心思想：传统NN是全连接的（所有输入直接连到所有隐藏层）。DAG-Net的网络结构由因果DAG决定——信息沿着因果路径流动。

标准NN:
  Input(47vars) → Hidden(全连接) → Output(presence)
  
DAG-Net:
  Layer 1: 只有DAG根节点(Elev, FlowAcc等) → Hidden₁
  Layer 2: Hidden₁ + 中间节点(BIO1-19) → Hidden₂  
  Layer 3: Hidden₂ + 叶节点(Soil, LC) → Output(presence)
  
  信息按照因果层级自上而下流动！
为什么好：

模型架构本身就编码了因果结构——这是CAST独有的
类似CISO用Transformer编码物种交互，CAST用DAG-Net编码变量间的因果交互
可以与CISO做方法论类比：CISO改造了输入（加biotic信息），CAST改造了架构（按因果结构布线）
难度：★★★★☆（需要自定义NN，但用torch或keras可实现）



步骤 1/8: 数据加载与因果层级分组...
  Layer 0 (Causal Roots):     6 vars [dem_range, dem_avg, slope_range, slope_avg, flow_length, flow_acc]
  Layer 1 (Causal Mediators): 19 vars
  Layer 2 (Causal Responses): 22 vars (12 LC + 10 Soil)
  Total: 47 vars

  ATE权重 (Top 5):
    lc_wavg_06: 1.0000
    lc_wavg_09: 0.1000
    soil_wavg_02: 0.1000
    lc_wavg_12: 0.1000
    flow_acc: 0.1000

步骤 2/8: 数据准备...
  训练: 4716 | 验证: 1179 | 测试: 1473

步骤 3/8: 定义DAG-Net (torch)...

步骤 4/8: 训练引擎...

步骤 5/8: 训练 DAG-Net Ensemble...

  === DAG-Net Run 1/5 ===
  Epoch  10: loss=0.7071 val_AUC=0.8809 lr=0.000997
  Epoch  20: loss=0.6658 val_AUC=0.8873 lr=0.000989
  Epoch  30: loss=0.6607 val_AUC=0.8907 lr=0.000976
  Epoch  40: loss=0.6393 val_AUC=0.8933 lr=0.000957
  Epoch  50: loss=0.6032 val_AUC=0.8951 lr=0.000934
  Epoch  60: loss=0.6013 val_AUC=0.8969 lr=0.000905
  Epoch  70: loss=0.5849 val_AUC=0.8970 lr=0.000873
  Epoch  80: loss=0.5835 val_AUC=0.9015 lr=0.000836
  Epoch  90: loss=0.5433 val_AUC=0.9011 lr=0.000796
  Epoch 100: loss=0.5173 val_AUC=0.8992 lr=0.000752
  Epoch 110: loss=0.5071 val_AUC=0.8938 lr=0.000706
  Early stop @ epoch 115 (best val_AUC=0.9015)
  Run 1: test_AUC=0.8898 val_AUC=0.9015

  === DAG-Net Run 2/5 ===
  Epoch  10: loss=0.6888 val_AUC=0.8809 lr=0.000997
  Epoch  20: loss=0.6613 val_AUC=0.8777 lr=0.000989
  Epoch  30: loss=0.6481 val_AUC=0.8827 lr=0.000976
  Epoch  40: loss=0.6267 val_AUC=0.8929 lr=0.000957
  Epoch  50: loss=0.6140 val_AUC=0.8930 lr=0.000934
  Epoch  60: loss=0.6044 val_AUC=0.8864 lr=0.000905
  Epoch  70: loss=0.5916 val_AUC=0.8888 lr=0.000873
  Epoch  80: loss=0.5627 val_AUC=0.8922 lr=0.000836
  Early stop @ epoch 82 (best val_AUC=0.8974)
  Run 2: test_AUC=0.8871 val_AUC=0.8974

  === DAG-Net Run 3/5 ===
  Epoch  10: loss=0.7069 val_AUC=0.8882 lr=0.000997
  Epoch  20: loss=0.6588 val_AUC=0.8835 lr=0.000989
  Epoch  30: loss=0.6553 val_AUC=0.8927 lr=0.000976
  Epoch  40: loss=0.6316 val_AUC=0.8813 lr=0.000957
  Epoch  50: loss=0.6214 val_AUC=0.8942 lr=0.000934
  Epoch  60: loss=0.6306 val_AUC=0.8883 lr=0.000905
  Epoch  70: loss=0.5936 val_AUC=0.8889 lr=0.000873
  Epoch  80: loss=0.5710 val_AUC=0.8933 lr=0.000836
  Epoch  90: loss=0.5648 val_AUC=0.8966 lr=0.000796
  Epoch 100: loss=0.5512 val_AUC=0.8967 lr=0.000752
  Epoch 110: loss=0.5498 val_AUC=0.8919 lr=0.000706
  Epoch 120: loss=0.5189 val_AUC=0.8923 lr=0.000658
  Epoch 130: loss=0.4997 val_AUC=0.8951 lr=0.000608
  Early stop @ epoch 136 (best val_AUC=0.8971)
  Run 3: test_AUC=0.8962 val_AUC=0.8971

  === DAG-Net Run 4/5 ===
  Epoch  10: loss=0.6789 val_AUC=0.8862 lr=0.000997
  Epoch  20: loss=0.6634 val_AUC=0.8928 lr=0.000989
  Epoch  30: loss=0.6424 val_AUC=0.8875 lr=0.000976
  Epoch  40: loss=0.6329 val_AUC=0.8906 lr=0.000957
  Epoch  50: loss=0.6024 val_AUC=0.8875 lr=0.000934
  Epoch  60: loss=0.6028 val_AUC=0.8905 lr=0.000905
  Epoch  70: loss=0.5936 val_AUC=0.8929 lr=0.000873
  Epoch  80: loss=0.5804 val_AUC=0.8961 lr=0.000836
  Epoch  90: loss=0.5559 val_AUC=0.8941 lr=0.000796
  Epoch 100: loss=0.5342 val_AUC=0.8966 lr=0.000752
  Early stop @ epoch 108 (best val_AUC=0.8984)
  Run 4: test_AUC=0.8928 val_AUC=0.8984

  === DAG-Net Run 5/5 ===
  Epoch  10: loss=0.7118 val_AUC=0.8799 lr=0.000997
  Epoch  20: loss=0.6753 val_AUC=0.8890 lr=0.000989
  Epoch  30: loss=0.6429 val_AUC=0.8876 lr=0.000976
  Epoch  40: loss=0.6167 val_AUC=0.8903 lr=0.000957
  Epoch  50: loss=0.6083 val_AUC=0.8905 lr=0.000934
  Epoch  60: loss=0.6009 val_AUC=0.8929 lr=0.000905
  Epoch  70: loss=0.5765 val_AUC=0.8919 lr=0.000873
  Epoch  80: loss=0.5795 val_AUC=0.8932 lr=0.000836
  Epoch  90: loss=0.5577 val_AUC=0.8947 lr=0.000796
  Epoch 100: loss=0.5505 val_AUC=0.8908 lr=0.000752
  Epoch 110: loss=0.5316 val_AUC=0.8907 lr=0.000706
  Early stop @ epoch 111 (best val_AUC=0.8982)
  Run 5: test_AUC=0.8917 val_AUC=0.8982

  Ensemble weights: 0.201 0.200 0.200 0.200 0.200 
  DAG-Net Ensemble: AUC=0.8953 TSS=0.6081
  DAG-Net Best:     AUC=0.8898 TSS=0.6041

步骤 6/8: Flat NN + Baselines...

  [Flat NN] Training...
  Epoch  10: loss=0.6357 val_AUC=0.8865
  Epoch  20: loss=0.6410 val_AUC=0.8895
  Epoch  30: loss=0.5956 val_AUC=0.8939
  Epoch  40: loss=0.5749 val_AUC=0.8897
  Epoch  50: loss=0.5290 val_AUC=0.8879
  Epoch  60: loss=0.5511 val_AUC=0.8920
  Epoch  70: loss=0.5320 val_AUC=0.8870
  Epoch  80: loss=0.5066 val_AUC=0.8920
  Epoch  90: loss=0.4910 val_AUC=0.8869
  Epoch 100: loss=0.4784 val_AUC=0.8856
  Epoch 110: loss=0.4744 val_AUC=0.8919
  Early stop @ epoch 110 (best=0.8974)
  Flat NN: AUC=0.8903 TSS=0.6194

  [CPF] Causal Prior Forest...
Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights.
  CPF: AUC=0.9520 TSS=0.7745

  [Maxent]...
  Maxent: AUC=0.9045 TSS=0.6352
  [RF]...
  RF: AUC=0.9566 TSS=0.7942

步骤 7/8: CAST Super-Ensemble...
  CAST Super: AUC=0.9301 TSS=0.6826

步骤 8/8: 结果汇总...

  ╔══════════════════════════════════════════════════════════════╗
  ║          Model Performance Comparison (Test Set)           ║
  ╠══════════════════════════════════════════════════════════════╣
  ║   Random Forest                  AUC=0.9566 TSS=0.7942 ║
  ║   Causal Prior Forest            AUC=0.9520 TSS=0.7745 ║
  ║ ★ CAST (DAG-Net + CPF)           AUC=0.9301 TSS=0.6826 ║
  ║   Maxent                         AUC=0.9045 TSS=0.6352 ║
  ║   DAG-Net Ensemble               AUC=0.8953 TSS=0.6081 ║
  ║   Flat NN                        AUC=0.8903 TSS=0.6194 ║
  ║   DAG-Net (Best Single)          AUC=0.8898 TSS=0.6041 ║
  ╚══════════════════════════════════════════════════════════════╝

  消融分析:
    DAG结构贡献:  ΔAUC = +0.0050 (DAG-Net Ensemble vs Flat NN)
    因果权重贡献: ΔAUC = -0.0046 (CPF vs RF)
  ✓ 图表已保存: figures/17_dagnet/performance_comparison.png

====================================================================== 
               DAG-Net 实验完成 (torch版)
======================================================================

  最优模型: Random Forest (AUC=0.9566, TSS=0.7942)

  输出文件:
    output/17_dagnet/comparison_table.csv
    figures/17_dagnet/performance_comparison.png

✓ 完成!