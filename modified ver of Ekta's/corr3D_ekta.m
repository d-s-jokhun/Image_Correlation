function value = corr3D_ekta(a,b)

T=size(a,1);

a_mean=mean(a);
b_mean=mean(b);
num=0;
denom_a=0;
denom_b=0;
for t=1:T
    if a(t,1)>-100000000 % so that blank spaces in the excel file will be ignored
    num=num+ dot ( (a(t,:)-a_mean) , (b(t,:)-b_mean) );
    denom_a=denom_a + dot ( (a(t,:)-a_mean) , (a(t,:)-a_mean) );
    denom_b=denom_b + dot ( (b(t,:)-b_mean) , (b(t,:)-b_mean) );
    end
end
denom=sqrt(denom_a)*sqrt(denom_b);
value = num/denom;