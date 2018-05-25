function Sxp=computaccuracy(trainsample,testsample)
test_tol=size(testsample,2);
train_tol=size(trainsample,2);
pre_label=zeros(1,test_tol);
h = waitbar(0,'Please wait...');
for i=1:test_tol
    xp=l1_ls(trainsample,testsample(:,i),0.01,0.01);
    per = i / test_tol;
    waitbar(per, h ,sprintf('%2.0f%%',per*100))
    Sxp(:,i)=xp;
end
close(h)

