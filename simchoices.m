clear

load ~/Dropbox/Projects/Ebbinghaus/code/testdataout.mat

S=50;

for s=1:S

X=cell2mat(data(s).X);
p=cell2mat(data(s).Z);

scale=1;
w(s)=.1 + s*0.02; %0.04
wx=.5;
sigma=1;

v=@(x,y) (x ./ (sigma + ((w(s)+wx)*x + w(s)*y)));

deu = p(1,:).*v(X(1,:),X(2,:)) - v(X(2,:),X(1,:));

%y=(scale*deu' + randn(1000,1) > 0);
y=(deu' + random('Logistic',0,scale,1000,1) > 0); %ensure logistic has variance scale^2

sum(y)

data(s).y=2-double(y);
data(s).W=mat2cell(repmat(s,1,1000),1,ones(1,1000));
end

save('~/Dropbox/Projects/Ebbinghaus/code/subjtestdataout2.mat','data')