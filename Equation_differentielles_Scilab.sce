fp=18;
fm=9;
mup=4;
mum=1;
K=5*10**3;
kpm=0.12;
kmp=0.043;
alpha=1.8;
steps=100000;
nema=5;
T=zeros(1,steps);
R=zeros(1,steps);
Rpm=zeros(1,steps);
Rp=zeros(1,steps);
Rm=zeros(1,steps);
P=zeros(nema,steps);
M=zeros(nema,steps);
p=zeros(nema,1);
m=zeros(nema,1);
Pinitial=100;
Minitial=0;
clf()
for k=1:nema do
p(k,1)=Pinitial;
m(k,1)=max(Minitial,1);
end
P(:,1)=p
M(:,1)=m
R(1)=sum((p+m).^2)/(sum(p)+sum(m)).^2
Rpm(1)=sum((p.*m))/(sum(p)*sum(m))
Rp(1)=sum((p).^2)/(sum(p)).^2
Rm(1)=sum((m).^2)/(sum(m)).^2
for j=2:steps do
        logi=1-(alpha*sum(p)+sum(m))/K
        tau=30*max(max(p),max(m))
        for k=1:nema do
            if rand() < fp*p(k)*logi/tau then
                if rand()<kpm then
                        m(k)=m(k)+1;
                else     
                        p(k)=p(k)+1;
                end
            end
            if rand() < fm*m(k)*logi/tau then
                if rand()<kmp then
                        p(k)=p(k)+1;
                else     
                        m(k)=m(k)+1;
                end
            end
            if rand() < mum*m(k)/tau then
                m(k)=m(k)-1;
            end
            if rand() < mup*p(k)/tau then
                p(k)=p(k)-1;
            end
        end
        T(j)=1/tau+T(j-1);
        for k=1:nema do
            P(k,j)=max(p(k),1);
            M(k,j)=max(m(k),1);
        end
        R(j)=sum((p+m).^2)/(sum(p)+sum(m)).^2
        Rpm(j)=sum((p.*m))/(sum(p)*sum(m))
        Rp(j)=sum((p).^2)/(sum(p)).^2
        Rm(j)=sum((m).^2)/(sum(m)).^2
    end
subplot(231)
plot2d("nl",T',[M',P'])
subplot(232)
plot2d(T,R)
subplot(233)
plot2d(T,Rpm)
subplot(234)
plot2d(T,Rm)
subplot(235)
plot2d(T,Rp)
