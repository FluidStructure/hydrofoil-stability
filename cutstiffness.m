c = 1;
for i = 1:30
    if (mod(i,6) == 2)||(mod(i,6) == 0)
        tmp1(:,c) = orgstiff(:,i);
        c = c+1;
    end
end

c = 1;
for i = 1:30
    if (mod(i,6) == 2)||(mod(i,6) == 0)
        tmp2(c,:) = tmp1(i,:);
        c = c+1;
    end
end

