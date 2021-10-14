%%%%%%%%%%四元数更新%%%%%%%%%%%
%% quaternion updata
%% 输入参数：T--仿真步长（秒）
%             Wibb--机体系陀螺仪输出   （单位：度/秒）
%             veloN--飞机运动速度——X东向、Y北向、Z天向（单位：米/秒）
%             posiN--经度、纬度、高度（单位：度、度、米）
%             Q--四元数
%% 输出参数：Q--四元数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q]=quaternion_updata(T,Wibb,Q)
Wnbb=Wibb*pi/180.0;  %单位：弧度/秒
%%%%%%%%%
WnbbA=Wnbb*T;
WnbbA0=sqrt(WnbbA(1,1)^2+WnbbA(2,1)^2+WnbbA(3,1)^2);
WnbbX=[0,          -WnbbA(1,1), -WnbbA(2,1), -WnbbA(3,1);
       WnbbA(1,1),  0,           WnbbA(3,1), -WnbbA(2,1);
       WnbbA(2,1), -WnbbA(3,1),   0,          WnbbA(1,1);
       WnbbA(3,1),  WnbbA(2,1),  -WnbbA(1,1),   0         ];
c_q=cos(WnbbA0/2);
if( WnbbA0<=1.0e-15 ) d_q=0.5; else  d_q=sin(WnbbA0/2)/WnbbA0;  end
Q=( c_q*eye(4)+d_q*WnbbX )*Q;%一般迭代算法
%%%%%%%%%%四元数规范化%%%%%%%%%
% tmp_Q=sqrt(Q(1,1)^2+Q(2,1)^2+Q(3,1)^2+Q(4,1)^2);
% for kc=1:4
%     Q(kc,1)=Q(kc,1)/tmp_Q;
% end
end

