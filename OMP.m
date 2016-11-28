function [A]=OMP(D,X,L)
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: 
%       D - the dictionary (its columns MUST be normalized).
%       X - the signals to represent
%       L - the max. number of coefficients for each signal.
% output arguments: 
%       A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);%n denotes the dimension of signals and P denotes the number of signals
[n,K]=size(D);%n denotes the dimension of dictionaries and P denotes the number of dictionaries
for k=1:1:P,%P denotes the number of signals
    a=[];
    x=X(:,k);%the k-th signal
    residual=x;
    indx=zeros(L,1);
    for j=1:1:L,
        proj=D'*residual;%find the predictor (column in A) most correlated with the residual
        [maxVal,pos]=max(abs(proj));%find the predictor (column in A) most correlated with the residual
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;%a denotes the optimal representation of ||x-D'a||^2_2, pinv denotes the Pseudoinverse
        residual=x-D(:,indx(1:j))*a;%update the residual
        if sum(residual.^2) < 1e-6
            break;
        end
    end;
    temp=zeros(K,1);
    temp(indx(1:j))=a;
    A(:,k)=sparse(temp);
end;
return;
