
Figure 1 :
  demo.m   
  DEBUG = 1&& SHOW = 1 为pipeline  输出为 saliency图  扩张凸包 原始凸包 3个cue prior_with_cue  种子点 
         
  DEBUG = 1&& SHOW = 0  输出为 saliency图   prior_with_cue  种子点

Figure 2:
 demo_illustion_convexhull.m  画凸包  扩张凸包 原始凸包
 gene_convex_curve.m

Figure 3：
  各种先验作为prior计算结果比较
  demo_illustion_compare_priors.m

Figure 4:
  使用test_sub_seeds.m
  
102 行   handles.database = 'MSRA5000'; % MSRA1000, Berkeley300
         用来选择需要的库 用的图要和这个库配套 
Figure 5:
demo_illustion_compare_seeds.m
改变  fseedOpts.fseed_num 用来选择前景点个数