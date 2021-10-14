function [ out_q ] = resample_AF(sigma_q,siran,particles,idd,mag,Fb )
current_q=particles;
g=9.7803698;
% [Xnext]=Artificial_fish_prey(try_number,X1,Y1,mag,Fb);
[random_q]=quaternion_sample(sigma_q);
sq=random_q;
m=1;
Cnb=[sq(2,m)^2+sq(1,m)^2-sq(4,m)^2-sq(3,m)^2, 2*(sq(2,m)*sq(3,m)+sq(1,m)*sq(4,m)), 2*(sq(2,m)*sq(4,m)-sq(1,m)*sq(3,m));
       2*(sq(2,m)*sq(3,m)-sq(1,m)*sq(4,m)), sq(3,m)^2-sq(4,m)^2+sq(1,m)^2-sq(2,m)^2,  2*(sq(3,m)*sq(4,m)+sq(1,m)*sq(2,m));
       2*(sq(2,m)*sq(4,m)+sq(1,m)*sq(3,m)), 2*(sq(3,m)*sq(4,m)-sq(1,m)*sq(2,m)), sq(4,m)^2-sq(3,m)^2-sq(2,m)^2+sq(1,m)^2]; 
siran_randomQ=mvnpdf([Cnb,zeros(3);zeros(3),Cnb]*[[0.48e-3;0.0;0.48e-3];0;0;-g],...
     [mag;Fb],1.*[2.5e-15*eye(3),zeros(3);zeros(3),8.37234225e-8*eye(3)]);
if((siran_randomQ/siran(idd))>1.03)
    current_q(:,idd)=random_q;
end
out_q=current_q;
end

