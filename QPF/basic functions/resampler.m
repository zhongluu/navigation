%%%%%%%%%%�ز���%%%%%%%%%%%
%%   resampler
%% ���������sin--��������--m��--m������
%            RegKG--��һ��������Ȩֵ--m��
%% ���������sout--�ز�������
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
    RegKG=ones(1,M)*(1/M);%��ʼ��Ȩֵ
end
    %% ����
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

