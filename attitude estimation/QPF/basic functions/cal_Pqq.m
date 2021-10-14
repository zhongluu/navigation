 %  Description:
 %   calculate quaternion second moment
 %  Reference:
 %  Crassidis, John L.Markley, F. Landis. Unscented Filtering for Spacecraft Attitude Estimation
 %  Oshman, Yaakov Carmi, Avishy.  Attitude Estimation from Vector Observations Using a Genetic-Algorithm-Embedded
 %                                     Quaternion Particle Filter
 %  Zhaihe ZHOU, Yulu ZHONG, Chuanwei ZENG, and Xiangrui TIAN. Attitude Estimation Using Parallel Quaternion Particle 
 %                                                                      Filter Based on New Quaternion Distribution
 % Declaration:
 %  Copyright(c) 2021-2025, by Yulu Zhong, Chuanwei Zeng, All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%% ���������sq--��Ԫ�����Ӽ���4��N�У�
%            RegKG--���Ӽ���Ӧ��Ȩֵ               
%% ���������Pqq--��Ԫ�����׾�
function [ Pqq ] = cal_Pqq( sq,RegKG )
%% 
[~,c]=size(sq);
Pqq=zeros(4,4);
for i=1:c
    
    Pqq=Pqq+RegKG(i)*sq(:,i)*sq(:,i)';

end

end

