%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    º½¼£·¢ÉúÆ÷
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025

function [t,atti,atti_rate]=traceset(t,T,atti,atti_rate)

if( t==0 )
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;
elseif(t<=5)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;
elseif(t<=6)
    atti_rate(1)=0.0;
    atti_rate(2)=2.0;
    atti_rate(3)=0.0;   
elseif(t<=9)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;  
elseif(t<=10)
    atti_rate(1)=0.0;
    atti_rate(2)=-2.0;
    atti_rate(3)=0.0;
elseif(t<=12)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=12.2)
    atti_rate(1)=0.0;
    atti_rate(2)=10.0;
    atti_rate(3)=10.0; 
elseif(t<=12.5)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;  
elseif(t<=12.7)
    atti_rate(1)=0.0;
    atti_rate(2)=-10.0;
    atti_rate(3)=-10.0; 
elseif(t<=13.0)
    atti_rate(1)=1.0*cos(2*pi*t);
    atti_rate(2)=1.0*cos(2*pi*t);
    atti_rate(3)=1.0*cos(2*pi*t); 
elseif(t<=15.0)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=15.5)
    atti_rate(1)=0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=20.0)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=20.5)
    atti_rate(1)=0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=21.5)
    atti_rate(1)=0.0;
    atti_rate(2)=2.0;
    atti_rate(3)=2.0; 
elseif(t<=23.5)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=24.5)
    atti_rate(1)=0.0;
    atti_rate(2)=-2.0;
    atti_rate(3)=-2.0; 
elseif(t<=29.5)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=30.0)
    atti_rate(1)=0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=34.5)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0;
%     atti_rate(3)=-2.0;
elseif(t<=35.0)
    atti_rate(1)=0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=36.0)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=39.0)
    atti_rate(1)=0.5*cos(2*pi*t);
    atti_rate(2)=0.5*cos(2*pi*t);
    atti_rate(3)=0.5*cos(2*pi*t);     
elseif(t<=40.0)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=41.0)
    atti_rate(1)=0.0;
    atti_rate(2)=-2.0;
    atti_rate(3)=0.0; 
elseif(t<=46.0)
    atti_rate(1)=0.0;
    atti_rate(2)=0.0;
    atti_rate(3)=0.0; 
elseif(t<=47.0)
    atti_rate(1)=0.0;
    atti_rate(2)=2.0;
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


