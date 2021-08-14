%%%%%%%%%%��Ԫ������%%%%%%%%%%%
%% quaternion updata
%% ���������T--���沽�����룩
%             Wibb--����ϵ���������   ����λ����/�룩
%             veloN--�ɻ��˶��ٶȡ���X����Y����Z���򣨵�λ����/�룩
%             posiN--���ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�
%             Q--��Ԫ��
%% ���������Q--��Ԫ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q]=quaternion_updata(T,Wibb,Q)
Wnbb=Wibb*pi/180.0;  %��λ������/��
%%%%%%%%%
WnbbA=Wnbb*T;
WnbbA0=sqrt(WnbbA(1,1)^2+WnbbA(2,1)^2+WnbbA(3,1)^2);
WnbbX=[0,          -WnbbA(1,1), -WnbbA(2,1), -WnbbA(3,1);
       WnbbA(1,1),  0,           WnbbA(3,1), -WnbbA(2,1);
       WnbbA(2,1), -WnbbA(3,1),   0,          WnbbA(1,1);
       WnbbA(3,1),  WnbbA(2,1),  -WnbbA(1,1),   0         ];
c_q=cos(WnbbA0/2);
if( WnbbA0<=1.0e-15 ) d_q=0.5; else  d_q=sin(WnbbA0/2)/WnbbA0;  end
Q=( c_q*eye(4)+d_q*WnbbX )*Q;%һ������㷨
%%%%%%%%%%��Ԫ���淶��%%%%%%%%%
% tmp_Q=sqrt(Q(1,1)^2+Q(2,1)^2+Q(3,1)^2+Q(4,1)^2);
% for kc=1:4
%     Q(kc,1)=Q(kc,1)/tmp_Q;
% end
end

