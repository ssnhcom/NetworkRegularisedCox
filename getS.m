% m: number of genes
% n: number of samples

function[S] = getNormalizedNetworkInformation(Data)

[m,n] = size(Data);
W = ones(m,m);
Data = Data';

Data = Data - ones(n,1)*mean(Data);
for i=1:m
    for j=i+1:m   
     c=sum(Data(:,i).*Data(:,j))/(sqrt(sum(Data(:,i).^2))*sqrt(sum(Data(:,j).^2)));
     W(i,j)=abs(c);
     W(j,i)=W(i,j);
    end
end

[~,IX]=sort(W,'descend');
[~,IXI]=sort(IX,'ascend');

clear c temp W IX 

W=zeros(m,m);
for i=1:m
    for j=1:m
        if i==j
        W(i,j) = 0;
        else
            W(i,j) = 1/(IXI(i,j))/(IXI(j,i));
        end
    end
end

Sum_R = sum(W, 2);
Sum_C = sum(W, 1);
S = zeros(m, m);
for i = 1 : m
    for j = 1 : m
            
        S(i, j) = W(i, j) / sqrt(Sum_R(i)) / sqrt(Sum_C(j));
        
    end
end

clear W;





