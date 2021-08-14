%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    º½¼£·¢ÉúÆ÷
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025

function [t,atti,atti_rate]=traceset4(t,T,atti,atti_rate)

if( t==0 )
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;
elseif(t<=0.1)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;
elseif(t<=5.100)
    atti_rate(1)=0.0;
    atti_rate(2)=15.0;
    atti_rate(3)=0.0;   
elseif(t<=9.000)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;  
elseif(t<=10.100)
    atti_rate(1)=0.0;
    atti_rate(2)=0;
    atti_rate(3)=0.0;
elseif(t<=12.000)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=13.000)
    atti_rate(1)=0;
    atti_rate(2)=1.0*cos(2*pi*t);
    atti_rate(3)=1.0*cos(2*pi*t); 
elseif(t<=15.100)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=20.100)
    atti_rate(1)=0.0;
    atti_rate(2)=-15.0;
    atti_rate(3)=20.0; 
elseif(t<=25.100)
    atti_rate(1)=0.0;
    atti_rate(2)=-15.0;
    atti_rate(3)=0.0; 
elseif(t<=30.000)
    atti_rate(1)=0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 

elseif(t<=35.000)
    atti_rate(1)=0.0;
    atti_rate(2)=15.0;
    atti_rate(3)=-20.0; 
elseif(t<=36.000)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=37.000)
    atti_rate(1)=0;
    atti_rate(2)=-1*cos(2*pi*t);
    atti_rate(3)=-1*cos(2*pi*t);     
elseif(t<=40.000)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=41.000)
    atti_rate(1)=0.0;
    atti_rate(2)=-30.0;
    atti_rate(3)=0.0; 
elseif(t<=46.000)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=47.000)
    atti_rate(1)=0.0;
    atti_rate(2)=30.0;
    atti_rate(3)=0.0; 
else
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
end

%    veloB(2,1)=veloB(2,1)+acceB(2,1)*T;    

   atti(1,1)=atti(1,1)+atti_rate(1,1)*T;
   atti(2,1)=atti(2,1)+atti_rate(2,1)*T;
   atti(3,1)=atti(3,1)+atti_rate(3,1)*T;
  
%   if(atti(1,1)>360)
%	atti(1,1)=atti(1,1)-fix(atti(1,1)/360)*360;
%   end


