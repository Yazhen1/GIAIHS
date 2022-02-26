function F = GIAIHS(M1,Pan, flag)
%Inputs:
%  Ml: low spatial resolution multispectral image;
%  Pan: high spatial reoslution panchormatic image;
%  beta: the parameter between 0 and 1, and used to balance 
%        the sharpness and smoothness, the smaller the smoothness;
%  flag: a constant that determine whether normalization or not, 
%        and 1 for normalization, 0 for not.
%Outputs:
%  F: fused image.
%--------------------------------------------------------------------------
% Implementation by:  Dr. Wang Yazhen  <15548300646@163.com>
% Please refer to the following paper(GIAIHS)
% Wang, Y.; Liu, G.; Zhang, R.; Liu, J. A Two-Stage Pansharpening Method for the Fusion of Remote-Sensing Images. Remote Sens. 2022, 14, 1121. https://doi.org/10.3390/rs14051121 
%
% At the same time, the code is modified from the IAIHS framework in the following literature
% Leung, Y.; Liu, J.M.; Zhang, J. An improved adaptive intensity hue saturation method for the fusion of remote sensing images.IEEE Geosci. Remote Sens. Lett.2013,11, 985–989
%
%--------------------------------------------------------------------------
if nargin<4
    flag = 1;
end

type = 0;
if nargin<3
    type = 1;
end
Mul= G1(M1,Pan);
[m,n,c] = size(Mul);
[p,q]   = size(Pan);

if m~=p && n~=q
    error('wrong input!');
end

Mul = double(Mul);
Pan = double(Pan);
%%%%%%%%%%
beta=0.15;
%%%%%%%%%%
if flag == 1
    %norm all bands of multspectral images to be in range 0 to 1.
    normcoeffs = ones(c,1);
    for i=1:c
        normcoeffs(i)=max(max(Mul(:,:,i)));
        Mul(:,:,i)=Mul(:,:,i)/normcoeffs(i);
    end
    %norm panchromatic image to be in range 0 to 1
    pancoeff = max(max(Pan));
    Pan=Pan/pancoeff;
end

%==========================================================================
alpha = findalph(Mul,Pan);
I = (alpha(1)*Mul(:,:,1)+alpha(2)*Mul(:,:,2)+alpha(3)*Mul(:,:,3)+alpha(4)*Mul(:,:,4));

%==========================================================================
%hist matching or not
Pan = (Pan-mean(Pan(:)))*std(I(:))/std(Pan(:)) + mean(I(:));
%end


lamda = 10^-9;
eps   = 10^-10;
Gpan  = expGdge(Pan,lamda,eps);

F = zeros(size(Mul));
for i=1:c
    Gmul = expGdge(Mul(:,:,i),lamda,eps);
    if type==1
        Temp = Mul(:,:,i); 
        beta = std(Temp(:))/(std(Temp(:))+std(I(:))); 
    end
    F(:,:,i)=Mul(:,:,i) + Mul(:,:,i)./mean(Mul,3).*(beta*Gpan+(1-beta)*Gmul).*(Pan-I);
end

if flag == 1
    %recover the multispectral images to their original range
    for i=1:c
        F(:,:,i)=F(:,:,i)*normcoeffs(i);
    end
end

%==========================================================================

function findalph = findalph(M, P)

[n, m, d] = size(M); %#ok

findalph = ones(d,1);

for i=1:d
    for j=1:d
        A(i,j) = sum(sum(M(:,:,i).*M(:,:,j)));
    end
    b(i,1) = sum(sum(P.*M(:,:,i)));
end

tau = 5;
iter = 150000;
gamma1 = 1/200000;
%gamma2 = 1;

inv = (eye(d) + 2*tau*gamma1*A)^(-1);

for i = 1:iter
    findalph = inv * (findalph+2*tau*max(-findalph,0)+2*tau*gamma1*b);
end

