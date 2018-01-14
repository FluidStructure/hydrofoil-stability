cntr = 1;
HALF_NPANELS = [10:10:120]
for j = HALF_NPANELS;
    j
    F = hydroIntegrals(0,0,0.5,'0015',j,1000,1,1);
    for i = 1:30
        Fmat(cntr,i) = F(i);
    end
    cntr = cntr + 1;
end

plot(HALF_NPANELS,Fmat)