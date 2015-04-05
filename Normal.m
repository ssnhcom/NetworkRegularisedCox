function[H] = Normal(H)

[n d] = size(H);

sum_t_H = sum(H,2);
sum_H = sum(H);


for i=1:n
   for j=1:d

     if(H(i,j)~=0)
         H(i,j) = H(i,j)./((sqrt(sum_t_H(i,1)))*(sqrt(sum_H(1,j))));
     end
   end
end
