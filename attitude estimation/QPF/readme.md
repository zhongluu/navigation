1、程序整体框架
---
主函数：quaternion_particle_filter.m
对比函数：SPF_of_atti.m、Copy_3_of_EKF.m、USQUE.m
对应论文SPF_of_atti.m-> Attitude Estimation from Vector Observations Using Genetic-Algorithm-Embedded Quaternion Particle Filter
Copy_3_of_EKF.m ->Kalman Filtering for Spacecraft Attitude Estimation
USQUE.m-> UNSCENTED FILTERING FOR SPACECRAFT ATTITUDE ESTIMATION
#
2、主要函数
---
particle_init：生成初始粒子群函数
Initalparticles：初始收敛函数->主要是改变过程噪声和量测噪声使误差逐步收敛并使用了KLD过程调节粒子数目
initalparticleswithoutKLD：为不使用KLD调节粒子数目的初始收敛函数
PQPF_alter:该函数为PQPF的可替换函数与PQPF效果相同
#
细节说明：
---
particle_init—使用观测矢量对初始化粒子见小论文3.1.  Quaternion Particle Initialization和论文” Attitude Estimation from Vector Observations Using Genetic-Algorithm-Embedded Quaternion Particle Filter”的附录部分。
	输入参数：b0—观测值 r0—参考值 M—粒子数
	输出参数：q—初始粒子集
initalparticleswithoutKLD --详见Particle Filtering for Sequential Spacecraft Attitude Estimation与普通粒子滤波不同点仅在于噪声变化，过程噪声变化的目的在于增加粒子多样性，量测噪声变化的目的在于防止粒子没有落入量测值范围内导致权值为零引起滤波发散。（KLD版本理论上来说鲁棒性要好一些但是时间消耗巨大，KLD引入的想法见论文”KLD-Sampling: Adaptive Particle Filters and Mobile Robot Localization”和” Adapting the Sample Size in Particle Filters Through KLD-Sampling”，KLD算法具体实现见论文” Adapting sample size in particle filters through KLD-resampling”—作者为T. Li, S. Sun和 T. Sattar）。
输入参数：t—时间 T—滤波周期 sq—粒子集 wibb—陀螺仪输出值 Fb—加速度计输出值 mag—磁力计输出值 RegKG—归一化权值 wibbAWR—陀螺漂移粒子集 noise_wARW—陀螺漂移白噪声大小
	输出参数：Pq—估计值 sq—粒子集 Neff—粒子有效性 RegKG—归一化权值 wibbAWR—陀螺漂移粒子集 M—粒子数
#
NOTE：We use BSD License
---
