%%%%%%%%%%重采样%%%%%%%%%%%
%%   resampler
%% 输入参数：sin--采样粒子--m列--m个粒子
%            RegKG--归一化后粒子权值--m列
%% 输出参数：sout--重采样粒子
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ sout,wARWout,RegKG ] = resampler(sin,wARWin,RegKG,threshold)
Neff=(1/sum(RegKG.^2));
[~,M]=size(sin);
sout=sin;
wARWout=wARWin;
if(Neff<threshold)
    u=rand(M,1);
    u=sort(u);
    l=cumsum(RegKG);
    i=1;
    for j=1:M
        while(i<=M)&&(u(i)<=l(j))
            sout(:,i)=sin(:,j);
            wARWout(:,i)=wARWin(:,j);
            i=i+1;
        end
    end
    RegKG=ones(1,M)*(1/M);%初始化权值
end
    %% 正则化
%     [ sout ] = regularize(sout,varPq );
% hopt=(4/(M*6))^(1/8);
% qvar=sqrt(varPq);
% qrnd=hopt*[qvar(1,1)*randn(1,M);qvar(2,2)*randn(1,M);qvar(3,3)*randn(1,M);qvar(4,4)*randn(1,M)];
% sout=sout+qrnd;
% qrnd=[ones(1,M);qrnd];
% for i=1:M
% sout(:,i)=qAntiMatrix(sout(:,i))*qrnd(:,i);
% end
end

