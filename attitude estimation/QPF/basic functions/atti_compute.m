%%%%%%%%%%%��̬����%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%% attitude compute
%% ���������Pq--��Ԫ��
%% ���������attiN--��̬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [attiN]=atti_compute(Q)
 %%%%%%%%%%��ȡ��̬����%%%%%%%%%
  Cbn=[Q(2,1)^2+Q(1,1)^2-Q(4,1)^2-Q(3,1)^2, 2*(Q(2,1)*Q(3,1)+Q(1,1)*Q(4,1)), 2*(Q(2,1)*Q(4,1)-Q(1,1)*Q(3,1));
       2*(Q(2,1)*Q(3,1)-Q(1,1)*Q(4,1)), Q(3,1)^2-Q(4,1)^2+Q(1,1)^2-Q(2,1)^2,  2*(Q(3,1)*Q(4,1)+Q(1,1)*Q(2,1));
       2*(Q(2,1)*Q(4,1)+Q(1,1)*Q(3,1)), 2*(Q(3,1)*Q(4,1)-Q(1,1)*Q(2,1)), Q(4,1)^2-Q(3,1)^2-Q(2,1)^2+Q(1,1)^2]; 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %����̬(���������������
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  attiN(1,1)=atan(-Cbn(1,3)/Cbn(3,3));
  attiN(2,1)=atan(Cbn(2,3)/sqrt(Cbn(2,1)*Cbn(2,1)+Cbn(2,2)*Cbn(2,2)));
  attiN(3,1)=atan(Cbn(2,1)/Cbn(2,2));
    %��λ������

  %�����ж�
  attiN(1,1)=attiN(1,1)*180.0/pi;
  attiN(2,1)=attiN(2,1)*180.0/pi;
  attiN(3,1)=attiN(3,1)*180.0/pi;
    % ��λ����

  if(Cbn(2,2)<0 ) 
   attiN(3,1)=180.0+attiN(3,1);
  else 
   if(Cbn(2,1)<0) attiN(3,1)=360.0+attiN(3,1); end
  end
    %����Ƕȣ���λ���ȣ�

  if(Cbn(3,3)<0)
   if(Cbn(1,3)>0) attiN(1,1)=180.0-attiN(1,1); end
   if(Cbn(1,3)<0) attiN(1,1)=-(180.0+attiN(1,1)); end
  end
    %����Ƕȣ���λ���ȣ�
end

