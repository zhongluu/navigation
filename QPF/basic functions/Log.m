function [ x_ba ] = Log( x,v )
[~,M]=size(x);x_ba=zeros(3,M);
for m=1:M
    alfa=acos(v'*x(:,m));
    x_=(x(:,m)-cos(alfa)*v)*alfa/sin(alfa);
    x_ba(:,m)=x_(2:4);
end
end

