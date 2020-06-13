%% testing figure 2
m = 100;
n = 1;
aa = rand(1,m)*2-1;
x = [aa,rand(1,n)+4];
y = [aa,rand(1,n)-5];
figure,plot(x,0.8*ones(length(x),1),'r.',y,0.2*ones(length(y),1),'b.');
set(gca,'YLim',[0,1]);
ux = mean(x);
uy = mean(y);
%%
t = 0;
xt = x+t;
uxt = mean(xt);
fx = exp(-0.5*(x-ux).^2)/sqrt(2*pi);
fxt = exp(-0.5*(xt-uxt).^2)/sqrt(2*pi);
fy = exp(-0.5*(y-uy).^2)/sqrt(2*pi);
figure,plot(1:length(x),fx,'r',1:length(x),fxt,'b',1:length(y),fy,'g');

