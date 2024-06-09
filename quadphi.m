function [C,m,s] = quadphi(A,ell_vector)
%EXPMHK  Matrix phik function.
%  2024/05/05 revised
%   [C,m,s] = ph_ll(A,ell_vector) computes the matrix phik of A 
%   using quadruple angle algorithm combined with Taylor series 
%   Inputs:
%       A: the input matrix
%       ell_vector = [i1,i2,i3]  the required indices of matrix function for computation.
%   Outputs: 
%      C:  the phi_{ell_vector} of matrix A.
%      m:  the approximation order used.
%      s:  the scaling parameter.
%   Reference:
%   [1]Dongping Li, Xue Yang, Xiuying Zhang. Efficient algorithm for oscillatory matrix functions. 2024.
%   Author: Dongping Li
   [m,s,pA,met]=select_m_s(A);
   ell=max(ell_vector);
   ell_leng=length(ell_vector);
   C=cell(1,ell_leng);
   if met==0
      C =tay_phi(pA,m,ell_vector);
      return;
   else
       pB=series_tays_phi(pA,m,s,ell_vector); 
       if s==0
        C=pB;
        return;
       end     
   end  
  P=[1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000, 51090942171709440000, 1124000727777607680000, 25852016738884978212864, 620448401733239409999872, 15511210043330986055303168, 403291461126605650322784256, 10888869450418351940239884288, 304888344611713836734530715648, 8841761993739700772720181510144, 265252859812191032188804700045312, 8222838654177922430198509928972288, 263130836933693517766352317727113216, 8683317618811885938715673895318323200, 295232799039604119555149671006000381952, 10333147966386144222209170348167175077888, 371993326789901177492420297158468206329856, 13763753091226343102992036262845720547033088, 523022617466601037913697377988137380787257344, 20397882081197441587828472941238084160318341120, 815915283247897683795548521301193790359984930816, 33452526613163802763987613764361857922667238129664, 1405006117752879788779635797590784832178972610527232, 60415263063373834074440829285578945930237590418489344, 2658271574788448529134213028096241889243150262529425408, 119622220865480188574992723157469373503186265858579103744, 5502622159812088456668950435842974564586819473162983440384, 258623241511168177673491006652997026552325199826237836492800, 12413915592536072528327568319343857274511609591659416151654400];   
  n=size(A,1);
  I=eye(n);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  S=cell(1,ell+1);
    for ii = 1:s-1 
        if ell==0
           S{1}=2*pB{1}*pB{1}-I; 
        else
           S{1}=2*pB{1}*pB{1}-I; 
           S{2}=pB{1}*pB{2};  
        for k=2:ell
                S{k+1}=pB{1}*pB{k+1}+pB{2}*pB{k};     
                for j=2:k
                    S{k+1}=S{k+1}+pB{j+1}/P(k-j+1);
                end
              S{k+1}=0.5^k*S{k+1}; 
        end
        end
        pB=S;
    end 
    %%%%%%%%%%%%%%%%%%%%%
    for i=1:ell_leng %最外层的输出
        h= ell_vector(i);
        if h==0  
           C{i}=2*pB{1}*pB{1}-I;
        elseif h==1 
           C{i}=pB{1}*pB{2};
        else
           C{i}=pB{1}*pB{h+1}+pB{2}*pB{h};     
           for j=2:h
               C{i}=C{i}+pB{j+1}/P(h-j+1);
           end
           C{i}=0.5^h*C{i};
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,s,pA,met]=select_m_s(A)
% [m,s,pA,met,im]=select_m_s_b(A)
% Given a square matrix A, this function obtains the order m and scaling
% parameter s for Taylor approximation (maximum order of m is 30)
% Optimal orders m:            2 4 6 9 12 16 20 25 30
% Corresponding matrix powers: 2 2 3 3  4  4  5  5  6
% Estimation of norms of matrix powers is used
%   Input:
%      A: the input matrix.
%
%   Outputs:
%      m: order or approximation to be used.
%      s: scaling parameter
%      pA: cell array with the powers of A, A^2,A^3,...,A^maxpow.
%      met:
%         met=0, Taylor approximation for Nilpotent matrices is used
%         met=1, Taylor approximation for not Nilpotent matrices is used
%
% Theta from forward and backward error bounds from [1] Table 2 
   mm=[2 4 6 9 12 16 20];
   P=[2 2 3 3 4  4 5];% m所需的中alpha_p(A)的可能p值
   Q=[2 2 3 3 4  4 5];% m所需的A的最高次幂  
   Theta=[5.16e-08   %m=1  forward bound
          4.31e-05  %m=2  forward bound
          1.32e-02  %m=4  forward bound
          1.92e-01  %m=6  forward bound
          1.75      %m=9  backward bound
          6.59      %m=12 backward bound
          21.09     %m=16 backward bound
          47.35];   %m=20 backward bound
%       Theta=[5.161913593731081e-08      %m=1  forward bound
%        4.307691256676447e-05      %m=2  forward bound
%        1.319680929892753e-02      %m=4  forward bound
%        1.895232414039165e-01      %m=6  forward bound
%        1.7985058769167590         %m=9  backward bound
%        6.752349007371135          %m=12 backward bound
%        9.971046342716772];        %m=16 backward bound
   s = 0;
   met=1;
   A=-A;
   pA{1}=A;d(1)=norm(A,1); 
   if d(1)==0,m=0;met=0; return;end %Null matrix
   if d(1)<=Theta(1),m=1;met=0;return;end %m=1  
   pA{2}=A*A;d(2)=norm(pA{2},1);  
   if d(2)==0,m=1;met=0;return;end %Nilpotent matrix
   a(1)=d(1);
   a(2)=max([d(2)^(1/2),(d(1)*d(2))^(1/3)]);
   eta=a(2);
   %m=2 
   if eta<=Theta(2)  
    m=2;
    return;
   end
   %m=4 
   if eta<=Theta(3)
    m=4;
    return;
   end   
   pA{3}=pA{2}*A;d(3)=norm(pA{3},1);
   if d(3)==0,m=2;met=0; return;end %Nilpotent matrix 
   a(2)=max([d(2)^(1/2),d(3)^(1/3)]);
   d(4)=min([d(1)*d(3),d(2)^2]);
   a(3)=max([d(3)^(1/3),d(4)^(1/4)]);   
   eta=min([a(2),a(3)]);
  % m = 6
  if eta<=Theta(4)
    m=6;
    return;
  end
% m = 9
  if eta<=Theta(5)
    m=9;
    return;
  end 
% m = 12
 pA{4}=pA{2}*pA{2};d(4)=norm(pA{4},1);
 if d(4)==0,m=3;met=0;return;end %Nilpotent matrix
 d(5)=min([d(1)*d(4),d(2)*d(3)]);
 a(3)=max([d(3)^(1/3),d(4)^(1/4)]);
 a(4)=max([d(4)^(1/4),d(5)^(1/5)]);
 eta=min([a(2),a(3),a(4)]);
if eta<=Theta(6)
    m=12;
    return;
end
% m = 16
if eta<=Theta(7)
    m=16;
    return;
end
% m = 20
 pA{5}=pA{4}*pA{1};d(5)=norm(pA{5},1); 
 if d(5)==0,m=4;met=0;return;end %Nilpotent matrix
 d(6)=min([d(1)*d(5),d(2)*d(4),d(3)*d(3)]);
 a(4)=max([d(4)^(1/4),d(5)^(1/5)]); 
 a(5)=max([d(5)^(1/5),d(6)^(1/6)]);
 eta=min([a(2),a(3),a(4),a(5)]);
 if eta<=Theta(8)
    m=20;
    return;
 end 
s= ceil(0.5*log2(eta/Theta(8)));
m=20;
for i=1:length(pA)
    pA{i}=pA{i}/4^(s*i);
end
end

function B=series_tays_phi(pA,m,s,ell_vector)
%   Computes the phi_{ik} of matrix A by a matrix Taylor method.
%   Inputs:
%     pA: cell array with the powers of A nedeed to compute Taylor
%          series, such that Ap{i} contains A^i, for i=2,3,...,q.
%      m:  order or approximation to be used.
%      s:  order of scaling
%      ell_vector: index of phi_ik needed to compute.
%   Outputs:
%      if s==0, output phi_{ik} of matrix A.
%      if s>=1, output phi_{k} of matrix A, k=0,1,...ell.
ell_leng=length(ell_vector);
ell=max(ell_vector);
n = size(pA{1},1);
pp =[1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000, 51090942171709440000, 1124000727777607680000, 25852016738884978212864, 620448401733239409999872, 15511210043330986055303168, 403291461126605650322784256, 10888869450418351940239884288, 304888344611713836734530715648, 8841761993739700772720181510144, 265252859812191032188804700045312, 8222838654177922430198509928972288, 263130836933693517766352317727113216, 8683317618811885938715673895318323200, 295232799039604119555149671006000381952, 10333147966386144222209170348167175077888, 371993326789901177492420297158468206329856, 13763753091226343102992036262845720547033088, 523022617466601037913697377988137380787257344, 20397882081197441587828472941238084160318341120, 815915283247897683795548521301193790359984930816, 33452526613163802763987613764361857922667238129664, 1405006117752879788779635797590784832178972610527232, 60415263063373834074440829285578945930237590418489344, 2658271574788448529134213028096241889243150262529425408, 119622220865480188574992723157469373503186265858579103744, 5502622159812088456668950435842974564586819473162983440384, 258623241511168177673491006652997026552325199826237836492800, 12413915592536072528327568319343857274511609591659416151654400];   
I = eye(n);
if s==0
 B=cell(1,ell_leng);
 switch m
    case 2 %q=2
        for j=1:ell_leng
        B{j}=zeros(n,n); 
        u=ell_vector(j)+1;
        p=pp(u:2:end);
        B{j}= B{j}+pA{2}/p(3) +pA{1}/p(2) + I/p(1);
        end     
	case 4 %q=2
        for j=1:ell_leng
        B{j}=zeros(n,n);
        u=ell_vector(j)+1;
        p=pp(u:2:end);
        B{j}=B{j}+ pA{2}/p(5)+ pA{1}/p(4) + I/p(3);
        B{j}=B{j}*pA{2}  + pA{1}/p(2)     + I/p(1);
        end
	case 6 %q=3
        for j=1:ell_leng
        B{j}=zeros(n,n);
        u=ell_vector(j)+1;
        p=pp(u:2:end);
        B{j}=B{j} + pA{3}/p(7)  + pA{2}/p(6) + pA{1}/p(5) + I/p(4);
        B{j}=B{j}*pA{3} + pA{2}/p(3)  + pA{1}/p(2)      + I/p(1);
        end
	case 9 %q=3
        for j=1:ell_leng
        B{j}=zeros(n,n);
        u=ell_vector(j)+1;
        p=pp(u:2:end);
        B{j}=B{j} + pA{3}/p(10)      + pA{2}/p(9) + pA{1}/p(8) + I/p(7);
        B{j}=B{j}*pA{3} + pA{2}/p(6)  + pA{1}/p(5) + I/p(4);
        B{j}=B{j}*pA{3} + pA{2}/p(3)  + pA{1}/p(2)  + I/p(1); 
        end
	case 12 %q=4
         for j=1:ell_leng
          B{j}=zeros(n,n);
          u=ell_vector(j)+1;
          p=pp(u:2:end);
          B{j}=B{j} + pA{4}/p(13)      + pA{3}/p(12) + pA{2}/p(11)+ pA{1}/p(10) + I/p(9);
          B{j}=B{j}*pA{4} + pA{3}/p(8)  + pA{2}/p(7)  + pA{1}/p(6) + I/p(5);
          B{j}=B{j}*pA{4} + pA{3}/p(4)  + pA{2}/p(3)  + pA{1}/p(2) + I/p(1);
        end
	 case 16 %q=4
          for j=1:ell_leng
          B{j}=zeros(n,n);
          u=ell_vector(j)+1;
          p=pp(u:2:end);
          B{j}=B{j} + pA{4}/p(17) + pA{3}/p(16) + pA{2}/p(15) + pA{1}/p(14) + I/p(13);
          B{j}=pA{4}*B{j} + pA{3}/p(12) + pA{2}/p(11) + pA{1}/p(10) + I/p(9);
          B{j}=pA{4}*B{j} + pA{3}/p(8)  + pA{2}/p(7)  + pA{1}/p(6)  + I/p(5);
          B{j}=pA{4}*B{j} + pA{3}/p(4)  + pA{2}/p(3)  + pA{1}/p(2)  + I/p(1);
          end
	 case 20 %q=5
          for j=1:ell_leng
          B{j}=zeros(n,n);
          u=ell_vector(j)+1;
          p=pp(u:2:end);
          B{j}=B{j}+ pA{5}/p(21)+ pA{4}/p(20) + pA{3}/p(19) + pA{2}/p(18) + pA{1}/p(17) + I/p(16);
          B{j}=pA{5}*B{j} + pA{4}/p(15) + pA{3}/p(14) + pA{2}/p(13) + pA{1}/p(12) + I/p(11);
          B{j}=pA{5}*B{j}+ pA{4}/p(10)  + pA{3}/p(9)  + pA{2}/p(8)  + pA{1}/p(7)  + I/p(6);
          B{j}=pA{5}*B{j}+ pA{4}/p(5)  + pA{3}/p(4)  + pA{2}/p(3)  + pA{1}/p(2)  + I/p(1);
         end
 end 
else
B=cell(1,ell+1);    
switch m
    case 2 %q=2
        for j=1:ell+1
        B{j}=zeros(n,n); 
        p=pp(j:2:end);
        B{j}= B{j}+pA{2}/p(3) +pA{1}/p(2) + I/p(1);
        end     
	case 4 %q=2
        for j=1:ell+1
        B{j}=zeros(n,n); 
        p=pp(j:2:end);
        B{j}=B{j}+ pA{2}/p(5)+ pA{1}/p(4) + I/p(3);
        B{j}=B{j}*pA{2}  + pA{1}/p(2)     + I/p(1);
        end
	case 6 %q=3
        for j=1:ell+1
        B{j}=zeros(n,n); 
        p=pp(j:2:end);
        B{j}=B{j} + pA{3}/p(7)  + pA{2}/p(6) + pA{1}/p(5) + I/p(4);
        B{j}=B{j}*pA{3} + pA{2}/p(3)  + pA{1}/p(2)      + I/p(1);
        end
	case 9 %q=3
        for j=1:ell+1
        B{j}=zeros(n,n); 
        p=pp(j:2:end);
        B{j}=B{j} + pA{3}/p(10)      + pA{2}/p(9) + pA{1}/p(8) + I/p(7);
        B{j}=B{j}*pA{3} + pA{2}/p(6)  + pA{1}/p(5) + I/p(4);
        B{j}=B{j}*pA{3} + pA{2}/p(3)  + pA{1}/p(2)  + I/p(1); 
        end
	case 12 %q=4
         for j=1:ell+1
          B{j}=zeros(n,n); 
          p=pp(j:2:end);
          B{j}=B{j} + pA{4}/p(13)      + pA{3}/p(12) + pA{2}/p(11)+ pA{1}/p(10) + I/p(9);
          B{j}=B{j}*pA{4} + pA{3}/p(8)  + pA{2}/p(7)  + pA{1}/p(6) + I/p(5);
          B{j}=B{j}*pA{4} + pA{3}/p(4)  + pA{2}/p(3)  + pA{1}/p(2) + I/p(1);
        end
	 case 16 %q=4
          for j=1:ell+1
          B{j}=zeros(n,n); 
          p=pp(j:2:end);
          B{j}=B{j} + pA{4}/p(17) + pA{3}/p(16) + pA{2}/p(15) + pA{1}/p(14) + I/p(13);
          B{j}=pA{4}*B{j} + pA{3}/p(12) + pA{2}/p(11) + pA{1}/p(10) + I/p(9);
          B{j}=pA{4}*B{j} + pA{3}/p(8)  + pA{2}/p(7)  + pA{1}/p(6)  + I/p(5);
          B{j}=pA{4}*B{j} + pA{3}/p(4)  + pA{2}/p(3)  + pA{1}/p(2)  + I/p(1);
          end
	 case 20 %q=5
          for j=1:ell+1
          B{j}=zeros(n,n); 
          p=pp(j:2:end);
          B{j}=B{j}+ pA{5}/p(21)+ pA{4}/p(20) + pA{3}/p(19) + pA{2}/p(18) + pA{1}/p(17) + I/p(16);
          B{j}=pA{5}*B{j} + pA{4}/p(15) + pA{3}/p(14) + pA{2}/p(13) + pA{1}/p(12) + I/p(11);
          B{j}=pA{5}*B{j}+ pA{4}/p(10)  + pA{3}/p(9)  + pA{2}/p(8)  + pA{1}/p(7)  + I/p(6);
          B{j}=pA{5}*B{j}+ pA{4}/p(5)  + pA{3}/p(4)  + pA{2}/p(3)  + pA{1}/p(2)  + I/p(1);
         end
end
end
end



function B= tay_phi(pA,m,ell_vector)
%[B,np]= tayr_exp(Ap,m)
%   Computes the exponential for nilpotent matrices  by Taylor
%   approximations.
%   Inputs:
%     pA: cell array with the powers of A nedeed to compute Taylor
%          series, such that pA{i} contains A^i, for i=2,3,...,q.
%      m:  order or approximation to be used.
%
%   Outputs:
%      B:  the exponential of matrix A.
%      np: matrix products carried out by the function.
pp =[1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000, 51090942171709440000, 1124000727777607680000, 25852016738884978212864, 620448401733239409999872, 15511210043330986055303168, 403291461126605650322784256, 10888869450418351940239884288, 304888344611713836734530715648, 8841761993739700772720181510144, 265252859812191032188804700045312, 8222838654177922430198509928972288, 263130836933693517766352317727113216, 8683317618811885938715673895318323200, 295232799039604119555149671006000381952, 10333147966386144222209170348167175077888, 371993326789901177492420297158468206329856, 13763753091226343102992036262845720547033088, 523022617466601037913697377988137380787257344, 20397882081197441587828472941238084160318341120, 815915283247897683795548521301193790359984930816, 33452526613163802763987613764361857922667238129664, 1405006117752879788779635797590784832178972610527232, 60415263063373834074440829285578945930237590418489344, 2658271574788448529134213028096241889243150262529425408, 119622220865480188574992723157469373503186265858579103744, 5502622159812088456668950435842974564586819473162983440384, 258623241511168177673491006652997026552325199826237836492800, 12413915592536072528327568319343857274511609591659416151654400];   
ell_leng=length(ell_vector);
B=cell(1,ell_leng);
n = size(pA{1},1); 
I = eye(n);
switch m 
    case 0
        for k=1:ell_leng
            j=ell_vector(k);
            B{k}=I/pp(j+1);
        end
	case 1
        for k=1:ell_leng
            j=ell_vector(k);
            B{k}= pA{1}/pp(j+3) + I/pp(j+1);
        end
	case 2 
        for k=1:ell_leng
            j=ell_vector(k);
            B{k}=pA{2}/pp(j+5) + pA{1}/pp(j+3)  + I/pp(j+1);
        end
    case 3  
        for k=1:ell_leng
         j=ell_vector(k);  
         B{k}=(pA{1}/pp(j+6)+I/pp(j+5))*pA{2} + pA{1}/pp(j+3) + I/pp(j+1);
        end
	case 4
        for k=1:ell_leng
        j=ell_vector(k);
        B{k}=pA{2}/pp(j+9) + pA{1}/pp(j+7)   + I/pp(j+5);
        B{k}=B{k}*pA{2}    + pA{1}/pp(j+3)   + I/pp(j+1);
        end
end
end

