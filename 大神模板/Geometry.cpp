const double eps=1e-8,pi=acos(-1.0);
//符号函数
int sgn(double x) {
    return (x>eps)-(x<-eps);
}
//平方函数
double sqr(double x) {
    return x*x;
}
/****************点****************/
struct node {
    double x,y;
    node() {}
    node(double a,double b):x(a),y(b) {}
//加
    node operator +(const node a)const {
        return node(x+a.x,y+a.y);
    }
//减
    node operator -(const node a)const {
        return node(x-a.x,y-a.y);
    }
//比例改变
    node operator /(const double a)const {
        return node(x/a,y/a);
    }
    void read() {
        scanf("%lf%lf",&x,&y);
    }
//两点距离
    double dis(node a) {
        return sqrt(sqr(x-a.x)+sqr(y-a.y));
    }
//向量旋转逆时针方向
    node tot(double a) {
        return node(x*cos(a)-y*sin(a),y*cos(a)+x*sin(a));
    }
//向量定长
    node trunc(double d) {
        double len=dis(node(0,0));
        return node(x*d/len,y*d/len);
    }
//重载小于
    bool operator <(const node a)const {
        if(sgn(x-a.x)==0) return sgn(y-a.y)<0;
        return sgn(x-a.x)<0;
    }
//重载等于
    bool operator ==(const node a)const {
        return !sgn(x-a.x)&&!sgn(y-a.y);
    }
};
//原点
node zero(0,0);
/****************圆****************/
struct circle {
    node cen;
    double r;
    void read() {
        cen.read();
        scanf("%lf",&r);
    }
//圆的面积
    double area() {
        return pi*sqr(r);
    }
} temp;
//判断a，b是否内含
bool inner(circle a,circle b) {
    if(sgn(a.r-b.r)>0) return 0;
    return sgn(a.cen.dis(b.cen)+a.r-b.r)<=0;
}
//按极角排序
bool cmp(const node a,const node b) {
    return atan2(a.y-temp.cen.y,a.x-temp.cen.x)<atan2(b.y-temp.cen.y,b.x-temp.cen.x);
}
//叉积
double mult(node p,node p1,node p2) {
    return (p1.x-p.x)*(p2.y-p.y)-(p1.y-p.y)*(p2.x-p.x);
}
//内积
double inner(node a,node b) {
    return a.x*b.x+a.y*b.y;
}
//向量夹角
double angle(node a,node b) {
    double temp=inner(a,b)/a.dis(zero)/b.dis(zero);
    if(sgn(temp-1)>=0) temp=1;
    if(sgn(temp+1)<=0) temp=-1;
    return acos(temp);
}
//弓形面积，按逆时针方向
double get(node a,node b) {
    double thera=angle(a-temp.cen,b-temp.cen);
    double s=sqr(temp.r)*thera/2.0;
    double s1=sqr(temp.r)*sin(thera)/2.0;
    if(sgn(mult(temp.cen,a,b))>0) return s-s1;
    return pi*sqr(temp.r)-s+s1;
}
//计算两圆交点（排除了相切的情况，不处理内含）
bool inter_or(circle c1,circle c2,node &p1,node &p2) {
    double len=c1.cen.dis(c2.cen);
    if(sgn(len-c1.r-c2.r)>=0) return 0;
    double s=(sqr(c1.r)-sqr(c2.r)+sqr(len))/len/2;
    double h=sqrt(sqr(c1.r)-sqr(s));
    node vec=c2.cen-c1.cen;
    node pp=c1.cen+vec.trunc(s);
    p1=pp+vec.tot(pi/2).trunc(h);
    p2=pp-vec.tot(pi/2).trunc(h);
    return 1;
}
//计算两圆交点（不处理内含）
bool inter_and(circle c1,circle c2,node &p1,node &p2) {
    double len=c1.cen.dis(c2.cen);
    if(sgn(len-c1.r-c2.r)>0) return 0;
    double s=(sqr(c1.r)-sqr(c2.r)+sqr(len))/len/2;
    double h=sqrt(sqr(c1.r)-sqr(s));
    node vec=c2.cen-c1.cen;
    node pp=c1.cen+vec.trunc(s);
    p1=pp+vec.tot(pi/2).trunc(h);
    p2=pp-vec.tot(pi/2).trunc(h);
    return 1;
}
/****************圆的集合****************/
struct H {
    int n;
    circle C[MAXN];
    void read() {
        for(int i=0; i<n; i++) C[i].read();
    }
//集合合并
    H operator +(const H a)const {
        H tmp;
        tmp.n=n+a.n;
        int i;
        for(i=0; i<n; i++) tmp.C[i]=C[i];
        for(i=n; i<tmp.n; i++) tmp.C[i]=a.C[i-n];
        return tmp;
    }
//圆的面积并去掉内含的圆
    void init_or() {
        bool mark[MAXN]= {0};
        int i,j,cnt=0;
        for(i=0; i<n; i++)
            for(j=0; j<n; j++) {
                if(i==j||mark[j]) continue;
                if(inner(C[i],C[j])) mark[i]=1;
            }
        for(i=0; i<n; i++) if(!mark[i]) C[cnt++]=C[i];
        n=cnt;
    }
//圆的面积交去掉内含的圆
    void init_and() {
        bool mark[MAXN]= {0};
        int i,j,cnt=0;
        for(i=0; i<n; i++)
            for(j=0; j<n; j++) {
                if(i==j||mark[j]) continue;
                if(inner(C[j],C[i])) mark[i]=1;
            }
        for(i=0; i<n; i++) if(!mark[i]) C[cnt++]=C[i];
        n=cnt;
    }
//判断圆弧是否被被其他圆盖住
    bool isvalid_or(node a,node b) {
        node c=(a+b)/2.0;
        node vec,p;
        int i;
        vec=a-b;
        p=temp.cen+vec.tot(pi/2).trunc(temp.r);
        for(i=0; i<n; i++) if(sgn(p.dis(C[i].cen)-C[i].r)<0) break;
        return i==n;
    }
//判断点是否被被其他圆盖住
    bool isvalid_and(node a) {
        int i;
        for(i=0; i<n; i++) if(sgn(a.dis(C[i].cen)-C[i].r)>0) return 0;
        return 1;
    }
//判断圆弧是否被被其他圆盖住
    bool isvalid_and(node a,node b) {
        node c=(a+b)/2.0;
        node vec,p;
        int i;
        vec=a-b;
        p=temp.cen+vec.tot(pi/2).trunc(temp.r);
        return isvalid_and(p);
    }
//圆的面积并
    double area_or() {
        init();
        int i,j,k;
        vector <node> S[MAXN];
        node a,b;
        for(i=0; i<n; i++)
            for(j=i+1; j<n; j++) {
                if(inter(C[i],C[j],a,b)) {
                    S[i].push_back(a);
                    S[i].push_back(b);
                    S[j].push_back(a);
                    S[j].push_back(b);
                }
            }
        double ans=0,tt=0;
        for(i=0; i<n; i++) {
            temp=C[i];
            if(S[i].size()==0) {
                ans+=C[i].area();
                continue;
            }
            sort(S[i].begin(),S[i].end(),cmp);
            S[i].resize(unique(S[i].begin(),S[i].end())-S[i].begin());
            if(S[i].front()==S[i].back()) S[i].pop_back();
            for(j=0; j<S[i].size(); j++) {
                k=(j+1)%S[i].size();
                if(isvalid_or(S[i][j],S[i][k])) {
                    ans+=get(S[i][j],S[i][k]);
                    tt+=mult(zero,S[i][j],S[i][k]);
                }
            }
        }
        return ans+fabs(tt/2.0);
    }
//圆的面积交
    double area_and(node &res) {
        init_and();
        int i,j,k;
        vector <node> S[MAXN];
        node a,b;
        for(i=0; i<n; i++)
            for(j=i+1; j<n; j++) {
                if(inter(C[i],C[j],a,b)) {
                    S[i].push_back(a);
                    S[i].push_back(b);
                    S[j].push_back(a);
                    S[j].push_back(b);
                } else return -1.0;
            }
        double ans=0,tt=0;
        for(i=0; i<n; i++) {
            temp=C[i];
            sort(S[i].begin(),S[i].end(),cmp);
            S[i].resize(unique(S[i].begin(),S[i].end())-S[i].begin());
            if(S[i].front()==S[i].back()) S[i].pop_back();
            if(S[i].size()==1) {
                if(isvalid_and(S[i][0])) res=S[i][0];
                continue;
            }
            for(j=0; j<S[i].size(); j++) {
                if(isvalid_and(S[i][j])) res=S[i][j];
                k=(j+1)%S[i].size();
                if(isvalid_and(S[i][j],S[i][k])) {
                    ans+=get(S[i][j],S[i][k]);
                    tt+=mult(zero,S[i][j],S[i][k]);
                }
            }
        }
        return ans+fabs(tt/2.0);
    }
}；
/****************多边形****************/
//极角排序
bool cmp(const node p1,const node p2) {
    if(sgn(mult(temp,p1,p2))==0) return temp.dis(p1)<temp.dis(p2);
    return sgn(mult(temp,p1,p2))>0;
}
struct poly {
    int n;
    node p[MAXN];
//多边形合并
    poly operator +(const poly a)const {
        poly tmp;
        tmp.n=n+a.n;
        for(i=0; i<n; i++) tmp.p[i]=p[i];
        for(i=n; i<tmp.n; i++) tmp.p[i]=a.p[i-n];
        return tmp;
    }
//判断点在多边形内
    bool in(node tmp) {
        int sum=0,i;
        for(i=0; i<n; i++) if(inter(tmp,I,p[i],p[(i+1)%n])) sum++;
        return sum&1;
    }
    void read() {
        for(int i=0; i<n; i++) p[i].read();
    }
//多边形面积
    double area() {
        double ans=0;
        for(int i=0; i<n; i++) ans+=mult(zero,p[i],p[(i+1)%n]);
        return fabs(ans/2);
    }
//凸包
    poly graham() {
        poly r;
        temp=node(inf,inf);
        int i,pos,top=2;
        for(i=0; i<n; i++) {
            if(temp.y>p[i].y) pos=i,temp=p[i];
            else if(sgn(temp.y-p[i].y)==0) {
                if(temp.x>p[i].x) pos=i,temp=p[i];
            }
        }
        p[pos]=p[0];
        p[0]=temp;
        sort(p+1,p+n,cmp);
        r.p[0]=p[0],r.p[1]=p[1];
        for(i=2; i<n; i++) {
            while(top>=2&&mult(r.p[top-2],r.p[top-1],p[i])<0) top--;
            if(top==1) {
                r.p[top++]=p[i];
                continue;
            }
            if(mult(r.p[top-2],r.p[top-1],p[i])>0) r.p[top++]=p[i];
            else if(r.p[top-2].dis(r.p[top-1])<r.p[top-2].dis(p[i])) r.p[top-1]=p[i];
        }
        r.n=top;
        return r;
    }
};
/****************小函数****************/
//点到直线距离
double linedis(node p,node p1,node p2) {
    return fabs(mult(p,p1,p2))/dis(p1,p2);
}
//点到线段距离
double segdis(node p,node p1,node p2) {
    node t=p;
    t.x+=p2.y-p1.y;
    t.y+=p1.x-p2.x;
    if(sgn(mult(p1,p,t)*mult(p2,p,t))>0) return dis(p,p1)<dis(p,p2)?dis(p,p1):dis(p,p2);
    return fabs(mult(p,p1,p2))/dis(p1,p2);
}
//求两直线交点
bool ins(node p1,node p2,node p3,node p4,node &tmp) {
    double a,b,c,a1,b1,c1;
    a=p1.y-p2.y;
    b=p2.x-p1.x;
    c=p2.y*p1.x-p1.y*p2.x;
    a1=p3.y-p4.y;
    b1=p4.x-p3.x;
    c1=p4.y*p3.x-p3.y*p4.x;
    if(sgn(a*b1-b*a1)==0) return 0;
    tmp=node((b*c1-b1*c)/(a*b1-b*a1),(a*c1-a1*c)/(a1*b-a*b1));
    return 1;
}
//两直线关系（所有情况）
//0一个交点，1无数个，2没有
int ins(node p1,node p2,node p3,node p4,node &tmp) {
    double a,b,c,a1,b1,c1;
    a=p1.y-p2.y;
    b=p2.x-p1.x;
    c=p2.y*p1.x-p1.y*p2.x;
    a1=p3.y-p4.y;
    b1=p4.x-p3.x;
    c1=p4.y*p3.x-p3.y*p4.x;
    if(sgn(a*b1-b*a1)==0) {
        if(sgn(a*c1-a1*c)==0&&sgn(b*c1-c*b1)==0) {
            return 1;
        }
        return 2;
    }
    tmp=node((b*c1-b1*c)/(a*b1-b*a1),(a*c1-a1*c)/(a1*b-a*b1));
    return 0;
}
//线段相交判断
bool check(node p1,node p2,node p3,node p4) {
    if(max(p1.x,p2.x)<min(p3.x,p4.x)) return 0;
    if(max(p1.y,p2.y)<min(p3.y,p4.y)) return 0;
    if(min(p1.x,p2.x)>max(p3.x,p4.x)) return 0;
    if(min(p1.y,p2.y)>max(p3.y,p4.y)) return 0;
    double aa=mult(p1,p2,p3);
    double bb=mult(p1,p2,p4);
    double cc=mult(p3,p4,p1);
    double dd=mult(p3,p4,p2);
    return sgn(aa*bb)<=0&&sgn(cc*dd)<=0;
}
//线段交点
bool ins(node p1,node p2,node p3,node p4,node &tmp) {
    double a=fabs(mult(p3,p1,p2));
    double b=fabs(mult(p4,p1,p2));
    if(!sgn(a)&&!sgn(b)) return 0;
    double tx=(p4.x*a+p3.x*b)/(a+b);
    double ty=(p4.y*a+p3.y*b)/(a+b);
    tmp=node(tx,ty);
    if(online(p1,p2,tmp)&&online(p3,p4,tmp)) return 1;
    return 0;
}
//判断点在线段上
bool online(node p1,node p2,node p) {
    return sgn(mult(p1,p2,p))==0&&sgn((p.x-p1.x)*(p.x-p2.x))<=0&&sgn((p.y-p1.y)*(p.y-p2.y))<=0;
}
//两圆面积交
double area(circle c1,circle c2) {
    double l=sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y));
    if(c1.r>=c2.r&&l<=c1.r-c2.r) return pi*c2.r*c2.r;
    if(c1.r<=c2.r&&l<=c2.r-c1.r) return pi*c1.r*c1.r;
    if(l>=c1.r+c2.r) return 0;
    double h1=(c1.r*c1.r-c2.r*c2.r+l*l)/(2*l);
    double len=sqrt(c1.r*c1.r-h1*h1);
    double h2=l-h1;
    double the1=acos(h1/c1.r);
    double the2=acos(h2/c2.r);
    return the1*c1.r*c1.r+the2*c2.r*c2.r-len*l;
}
//三角形外接圆
node outercircle(node tri[]) {
    double a,b,c,xA,xB,xC,c1,c2,yA,yB,yC;
    a=tri[0].dis(tri[1]);
    b=tri[1].dis(tri[2]);
    c=tri[2].dis(tri[0]);
    cr=a*b*c/fabs(mult(tri[0],tri[1],tri[2]))/2;
    xA = tri[0].x;
    yA = tri[0].y;
    xB = tri[1].x;
    yB = tri[1].y;
    xC = tri[2].x;
    yC = tri[2].y;
    c1 = (xA * xA + yA * yA - xB * xB - yB * yB) / 2;
    c2 = (xA * xA + yA * yA - xC * xC - yC * yC) / 2;

    circle.x = (c1 * (yA - yC) - c2 * (yA - yB)) /
               ((xA - xB) * (yA - yC) - (xA - xC) * (yA - yB));
    circle.y = (c1 * (xA - xC) - c2 * (xA - xB)) /
               ((yA - yB) * (xA - xC) - (yA - yC) * (xA - xB));
}
//两点求圆
bool get_center(node p1,node p2,node A[]) {
    double len=p1.dis(p2);
    node tmp=(p1+p2)/2;
    if(sgn(len-2*r)>0) return 0;
    double s=sqrt(sqr(r)-sqr(len/2));
    node vec=(p2-p1).rot(pi/2).trunc(s);
    A[0]=tmp+vec,A[1]=tmp-vec;
    return 1;
}
//一线一点求圆
bool get_center(node p,node p1,node p2,node A[]) {
    double len=fabs(mult(p,p1,p2)/p1.dis(p2));
    if(sgn(len-2*r)>0||sgn(len)==0) return 0;
    double s=sqrt(sqr(r)-sqr(r-len));
    node vec=(p2-p1).trunc(s);
    node V=vec.rot(pi/2);
    if(sgn(dot(V,p-p1))<0) V=zero-V;
    V=V.trunc(r-len);
    A[0]=p+vec+V;
    A[1]=p-vec+V;
    return 1;
}
//两线求圆
bool get_center(node p1,node p2,node p3,node p4,node A[]) {
    node tmp;
    if(!ins(p1,p2,p3,p4,tmp)) return 0;
    node v1=(p2-p1).rot(pi/2).trunc(r),v2=(p4-p3).rot(pi/2).trunc(r);
    int i;
    for(i=0; i<=3; i++) ins(p1+v1*(double)dx[i],p2+v1*(double)dx[i],p3+v2*(double)dy[i],p4+v2*(double)dy[i],A[i]);
    return 1;
}
//圆与直线交点
int cir_seg(node cir,double r,node a,node b,node &p1,node &p2) {
    node vec=b-a;
    node tmp;
    ins(a,b,cir,cir+vec.tot(pi/2),tmp);
    if(r-cir.dis(tmp)<0) return 0;
    double len=sqrt(sqr(r)-sqr(cir.dis(tmp)));
    p1=tmp-vec.trunc(len);
    p2=tmp+vec.trunc(len);
    int ans=0;
    if(online(a,b,p1)) ans++;
    if(online(a,b,p2)) ans++;
    if(online(a,b,p2)&&!online(a,b,p1)) tmp=p1,p1=p2,p2=tmp;
    return ans;
}
//最小圆覆盖
double cr;
int n;
void find(int state,node tri[]) {
    if(state==0) {
        circle=node(0,0);
        cr=-1;
    }
    if(state==1) {
        circle=tri[0];
        cr=0;
    }
    if(state==2) {
        circle=node((tri[0].x+tri[1].x)/2,(tri[0].y+tri[1].y)/2);
        cr=tri[0].dis(tri[1])/2;
    }
    if(state==3) outercircle(tri);
}
void mincircle(int cur,int state,node tri[]) {
    find(state,tri);
    if(state==3) return ;
    for(int i=1; i<=cur; i++) {
        if(dis(rec[i],circle)>cr) {
            tri[state]=rec[i];
            mincircle(i-1,state+1,tri);
        }
    }
}
void solve() {
    node tri[4];
    mincircle(n,0,tri);
    printf("%.2lf %.2lf %.2lf\n",circle.x,circle.y,cr);
}
int main() {
    while(scanf("%d",&n),n) {
        for(int i=1; i<=n; i++)
            scanf("%lf%lf",&rec[i].x,&rec[i].y);
        solve();
    }
}
//半平面交(p1,p2顺时针方向)
void cut(node p1,node p2) {
    int i,j;
    poly temp;
    node tmp;
    double a,b;
    temp.n=0;
    for(i=0; i<rec.n; i++) {
        j=(i+1)%rec.n;
        a=mult(p1,p2,rec.p[i]);
        b=mult(p1,p2,rec.p[j]);
        if(sgn(a)<=0) temp.p[temp.n++]=rec.p[i];
        if(sgn(a*b)<0) {
            ins(p1,p2,rec.p[i],rec.p[j],tmp);
            temp.p[temp.n++]=tmp;
        }
    }
    rec=temp;
}
//半平面交求核O(nlogn)
struct node {
    double x,y;
    node() {}
    node(double a,double b):x(a),y(b) {}
} rec[20010];
struct Line {
    node l,r;
} line[20010];
double at[20010];
int dq[20010],ord[20010];
int n;
bool judge(int a,int b,int c) {
    node temp;
    ins(line[a].l,line[a].r,line[b].l,line[b].r,temp);
    return mult(line[c].l,line[c].r,temp)<0;
}
bool cmp(const int a,const int b) {
    if(fabs(at[a]-at[b])<eps) return mult(line[b].l,line[b].r,line[a].l)>0;
    return at[a]<at[b];
}
void solve() {
    int i,j,bot,top;
    double x2,y1,x1,y2,temp,ans;
    node tmp;
    for(i=0; i<n; i++) {
        scanf("%lf%lf%lf%lf",&x1,&y1,&x2,&y2);
        line[i].l=node(x1,y1);
        line[i].r=node(x2,y2);
    }
    line[n].l=node(0,0);
    line[n++].r=node(inf,0);
    line[n].l=node(inf,0);
    line[n++].r=node(inf,inf);
    line[n].l=node(inf,inf);
    line[n++].r=node(0,inf);
    line[n].l=node(0,inf);
    line[n++].r=node(0,0);
    for(i=0; i<n; i++) {
        at[i]=atan2(line[i].r.y-line[i].l.y,line[i].r.x-line[i].l.x);
        ord[i]=i;
    }
    sort(ord,ord+n,cmp);
    temp=at[ord[0]];
    for(i=1,j=1; i<n; i++) {
        if(sgn(temp-at[ord[i]])==0) continue;
        ord[j++]=ord[i];
        temp=at[ord[i]];
    }
    n=j;
    dq[bot=1]=ord[0];
    dq[top=2]=ord[1];
    for(i=2; i<n; i++) {
        while(bot<top&&judge(dq[top-1],dq[top],ord[i])) top--;
        while(bot<top&&judge(dq[bot+1],dq[bot],ord[i])) bot++;
        dq[++top]=ord[i];
    }
    while(bot<top&&judge(dq[top-1],dq[top],dq[bot])) top--;
    while(bot<top&&judge(dq[bot+1],dq[bot],dq[top])) bot++;
    if(top-bot<2) {
        printf("0.0\n");
        return ;
    }
    dq[--bot]=dq[top];
    j=0;
    for(i=bot; i<top; i++) {
        ins(line[dq[i]].l,line[dq[i]].r,line[dq[i+1]].l,line[dq[i+1]].r,tmp);
        rec[j++]=tmp;
    }
    ans=0;
    for(i=0; i<j; i++)
        ans+=mult(node(0,0),rec[i],rec[(i+1)%j]);
    printf("%.1lf\n",fabs(ans/2));
}
int main() {
    while(~scanf("%d",&n)) solve();
}

//三角形与圆交 by ZaakDov
struct stct0 {
    double x,y;
};
typedef struct stct0 point;
double outer(point a,point b,point c) {
    return (a.x-c.x)*(b.y-c.y)-(a.y-c.y)*(b.x-c.x);
}
double inner(point a,point b,point c) {
    return (a.x-c.x)*(b.x-c.x)+(a.y-c.y)*(b.y-c.y);
}
double dist(point a,point b) {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
double calcarea(point a,point b,point c,double r) {
    double A,B,C,x,y,tS;
    A=dist(b,c);
    B=dist(a,c);
    C=dist(b,a);
    if(A<r&&B<r)
        return outer(a,b,c)/2;
    else if(A<r&&B>=r) {
        x=(inner(a,c,b)+sqrt(r*r*C*C-outer(a,c,b)*outer(a,c,b)))/C;
        tS=outer(a,b,c)/2;
        return asin(tS*(1-x/C)*2/r/B)*r*r/2+tS*x/C;
    } else if(A>=r&&B<r) {
        y=(inner(b,c,a)+sqrt(r*r*C*C-outer(b,c,a)*outer(b,c,a)))/C;
        tS=outer(a,b,c)/2;
        return asin(tS*(1-y/C)*2/r/A)*r*r/2+tS*y/C;
    } else if(fabs(outer(a,b,c))>=r*C||inner(b,c,a)<=0||inner(a,c,b)<=0)
        if(inner(a,b,c)<0)
            if(outer(a,b,c)<0)
                return (-acos(-1)-asin(outer(a,b,c)/A/B))*r*r/2;
            else
                return (acos(-1)-asin(outer(a,b,c)/A/B))*r*r/2;
        else
            return asin(outer(a,b,c)/A/B)*r*r/2;
    else {
        x=(inner(a,c,b)+sqrt(r*r*C*C-outer(a,c,b)*outer(a,c,b)))/C;
        y=(inner(b,c,a)+sqrt(r*r*C*C-outer(b,c,a)*outer(b,c,a)))/C;
        tS=outer(a,b,c)/2;
        return (asin(tS*(1-x/C)*2/r/B)+asin(tS*(1-y/C)*2/r/A))*r*r/2+tS*((y+x)/C-1);
    }
}
main() {
    point p[3],cen;
    double r,total;
    int i;
    while(~scanf("%lf%lf%lf",&cen.x,&cen.y,&r)) {
        for(i=0; i<=2; i++)
            scanf("%lf%lf",&p[i].x,&p[i].y);
        total=0;
        for(i=0; i<=2; i++)
            total+=calcarea(p[i],p[(i+1)%3],cen,r);
        printf("%.2lf\n",fabs(total));
    }
    return 0;
}
//三角形与圆交 by laprovence
double sum(node cir,double r,node a,node b) {
    double ra=cir.dis(a),rb=cir.dis(b);
    double ans=0,the,len;
    node ta,tb,tc,td,tt;
    double oper;
    if(mult(cir,a,b)==0) return 0;
    if(mult(cir,a,b)<0) oper=-1;
    else oper=1;
    if(ra<=r&&rb<=r) ans=fabs(mult(cir,a,b)/2);
    else if((ra<r&&rb>r)||(rb<r&&ra>r)) {
        if(ra<rb) cir_seg(cir,r,cir,b,ta,tt);
        else cir_seg(cir,r,cir,a,ta,tt);
        cir_seg(cir,r,a,b,tc,tt);
        if(ra<rb) ans+=fabs(mult(cir,a,tc)/2);
        else ans+=fabs(mult(cir,b,tc)/2);
        len=ta.dis(tc);
        the=acos((r*r*2-len*len)/r/r/2);
        ans+=sqr(r)*the/2;
    } else {
        if(cir_seg(cir,r,a,b,ta,tb)==2) {
            ans+=fabs(mult(cir,ta,tb)/2);
            cir_seg(cir,r,cir,a,tc,tt);
            len=ta.dis(tc);
            the=acos((r*r*2-len*len)/r/r/2);
            ans+=sqr(r)*the/2;
            cir_seg(cir,r,cir,b,td,tt);
            len=tb.dis(td);
            the=acos((r*r*2-len*len)/r/r/2);
            ans+=sqr(r)*the/2;
        } else {
            cir_seg(cir,r,cir,a,ta,tt);
            cir_seg(cir,r,cir,b,tb,tt);
            len=ta.dis(tb);
            the=acos((r*r*2-len*len)/r/r/2);
            ans+=sqr(r)*the/2;
        }
    }
    return ans*oper;
}
int main() {
    int i;
    double r,ans;
    while(~scanf("%lf%lf",&P[0].x,&P[0].y)) {
        for(i=1; i<=3; i++) P[i].read();
        scanf("%lf",&r);
        ans=0;
        for(i=0; i<=2; i++) ans+=sum(P[3],r,P[i],P[(i+1)%3]);
        printf("%.2lf\n",fabs(ans)+eps);
    }
    return 0;
}

/****************三维计算几何模板****************/
//空间中的点
struct node {
    double x,y,z;
    node() {}
    node(double a,double b,double c):x(a),y(b),z(c) {}
    void read() {
        scanf("%lf%lf%lf",&x,&y,&z);
    }
    node operator -(node a) {
        return node(x-a.x,y-a.y,z-a.z);
    }
    node operator *(node a) {
        return node(y*a.z-a.y*z,a.x*z-x*a.z,x*a.y-a.x*y);
    }
    double operator ^(node a) {
        return x*a.x+y*a.y+z*a.z;
    }
    double dis(node a) {
        return sqrt(sqr(x-a.x)+sqr(y-a.y)+sqr(z-a.z));
    }
} zero(0,0,0);
//向量叉积
node xmult(node u,node v) {
    return node(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x);
}
//计算向量
node sub(node u,node v) {
    return node(v.x-u.x,v.y-u.y,v.z-u.z);
}
//计算平面法向量
node pvec(node p1,node p2,node p3) {
    return xmult(sub(p1,p2),sub(p1,p3));
}
//计算向量内积
double dot(node p1,node p2) {
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}
//计算两点距离
double dis(node p1,node p2) {
    return sqrt(sqr(p1.x-p2.x)+sqr(p1.y-p2.y)+sqr(p1.z-p2.z));
}
//计算直线与平面交点
node inter(node p1,node p2,node p3,node u,node v) {
    node ret=pvec(p1,p2,p3),zero(0,0,0);
    double t=dot(ret,sub(u,p1))/dot(ret,sub(u,v));
    node vec=sub(u,v);
    ret.x=u.x+vec.x*t;
    ret.y=u.y+vec.y*t;
    ret.z=u.z+vec.z*t;
    return ret;
}
//点关于任意轴的旋转矩阵（不包括竖直的情况）（x，y，z，1）右乘该矩阵
void get_matri(point start,point end,double thera,double matri[][5]) {
    double tmp[5][5];
    double a[4],b[4],c[4],D[4],d[4],ds;
    double x1=start.x,x2=end.x,y1=start.y,y2=end.y,z1=start.z,z2=end.z;
    ds=sqrt(sqr(end.x-start.x)+sqr(end.y-start.y)+sqr(end.z-start.z));
    a[1]=(end.x-start.x)/ds;
    b[1]=(end.y-start.y)/ds;
    c[1]=(end.z-start.z)/ds;
    a[2]=(-b[1])/sqrt(sqr(a[1])+sqr(b[1]));
    b[2]=a[1]/sqrt(sqr(a[1])+sqr(b[1]));
    c[2]=0;
    a[3]=-b[2]*c[1],b[3]=a[2]*c[1],c[3]=sqrt(sqr(a[1])+sqr(b[1]));
    d[1]=-(a[1]*x1+b[1]*y1+c[1]*z1);
    d[2]=-(a[2]*x1+b[2]*y1+0*z2);
    d[3]=-(a[3]*x1+b[3]*y1+c[3]*z1);
    D[1]=-(a[1]*d[1]+a[2]*d[2]+a[3]*d[3]);
    D[2]=-(b[1]*d[1]+b[2]*d[2]+b[3]*d[3]);
    D[3]=-(c[1]*d[1] +c[3]*d[3]);
    int i,j,k;
    for(i=1; i<=4; i++) {
        for(j=1; j<=4; j++) {
            if(j==4) {
                if(i<4) matri[i][j]=0;
                else matri[i][j]=1;
                continue;
            }
            if(i==1) matri[i][j]=a[j];
            if(i==2) matri[i][j]=b[j];
            if(i==3) matri[i][j]=c[j];
            if(i==4) matri[i][j]=d[j];
        }
    }
    for(i=1; i<=4; i++)
        for(j=1; j<=4; j++)
            tmp[i][j]=0;
    tmp[1][1]=1,tmp[4][4]=1,tmp[2][2]=cos(thera);
    tmp[2][3]=sin(thera),tmp[3][2]=-sin(thera),tmp[3][3]=cos(thera);
    mul(matri,tmp);
    for(i=1; i<=4; i++)
        for(j=1; j<=4; j++)
            tmp[i][j]=0;
    for(i=1; i<=3; i++)
        for(j=1; j<=3; j++) {
            if(i==1) tmp[j][i]=a[j];
            if(i==2) tmp[j][i]=b[j];
            if(i==3) tmp[j][i]=c[j];
        }
    for(i=1; i<=3; i++) tmp[4][i]=D[i];
    tmp[4][4]=1;
    mul(matri,tmp);
}
//三维凸包
int to[MAXN][MAXN];
struct ThreeD {
    int n;
    node p[MAXN];
    struct Face {
        int a,b,c;
        bool ok;
    } F[MAXN*8];
    int cnt;
    void read() {
        for(int i=0; i<n; i++) p[i].read();
    }
    double volume(node p1,node p2,node p3,node p4) {
        return ((p2-p1)*(p3-p1))^(p4-p1);
    }
    double get(node p1,Face f) {
        return volume(p[f.a],p[f.b],p[f.c],p1);
    }
    void deal(int pos,int a,int b) {
        int f=to[a][b];
        Face tmp;
        if(F[f].ok) {
            if(sgn(get(p[pos],F[f]))>0) dfs(pos,f);
            else {
                to[b][a]=to[a][pos]=to[pos][b]=cnt;
                tmp.a=b,tmp.b=a,tmp.c=pos,tmp.ok=1;
                F[cnt++]=tmp;
            }
        }
    }
    void dfs(int pos,int f) {
        F[f].ok=0;
        deal(pos,F[f].b,F[f].a);
        deal(pos,F[f].c,F[f].b);
        deal(pos,F[f].a,F[f].c);
    }
    void run() {
        cnt=0;
        if(n<=3) return ;
        int i;
        for(i=1; i<n; i++)
            if(sgn(p[0].dis(p[i]))>0) {
                swap(p[1],p[i]);
                break;
            }
        if(i==n) return ;
        for(i=2; i<n; i++) {
            if(sgn(zero.dis((p[0]-p[1])*(p[i]-p[0])))>0) {
                swap(p[2],p[i]);
                break;
            }
        }
        if(i==n) return ;
        for(i=3; i<n; i++) {
            if(fabs(volume(p[0],p[1],p[2],p[i]))>eps) {
                swap(p[3],p[i]);
                break;
            }
        }
        if(i==n) return ;
        Face tmp;
        for(i=0; i<=3; i++) {
            tmp.a=(i+1)%4,tmp.b=(i+2)%4,tmp.c=(i+3)%4,tmp.ok=1;
            if(sgn(get(p[i],tmp))>0) swap(tmp.b,tmp.c);
            to[tmp.a][tmp.b]=to[tmp.b][tmp.c]=to[tmp.c][tmp.a]=cnt;
            F[cnt++]=tmp;
        }
        int j;
        for(i=4; i<n; i++)
            for(j=0; j<cnt; j++) {
                if(F[j].ok&&sgn(get(p[i],F[j]))>0) {
                    dfs(i,j);
                    break;
                }
            }
        int T=cnt;
        cnt=0;
        for(i=0; i<T; i++) if(F[i].ok) F[cnt++]=F[i];
    }
    bool same(Face A,Face B) {
        return sgn(get(p[A.a],B))==0&&
               sgn(get(p[A.b],B))==0&&
               sgn(get(p[A.c],B))==0;
    }
    double area() {
        int i;
        double ans=0;
        for(i=0; i<cnt; i++) ans+=zero.dis((p[F[i].b]-p[F[i].a])*(p[F[i].c]-p[F[i].a]))/2.0;
        return ans;
    }
    int FaceNum() {
        int i,j,ans=0;
        for(i=0; i<cnt; i++) {
            for(j=i+1; j<cnt; j++) {
                if(same(F[i],F[j])) break;
            }
            if(j==cnt) ans++;
        }
        return ans;
    }
} H;
int main() {
    while(~scanf("%d",&H.n)) {
        H.read();
        H.run();
        printf("%.3lf\n",H.area()+eps);
    }
    return 0;
}

//任意多边形的面积交（两个）
#include<iostream>
#include<stdio.h>
#include<math.h>
using namespace std;
const double eps=1e-8;
int sgn(double x) {
    return (x>eps)-(x<-eps);
}
struct node {
    double x,y;
    node() {}
    node(double a,double b):x(a),y(b) {}
    void read() {
        scanf("%lf%lf",&x,&y);
    }
};
double mult(node p,node p1,node p2) {
    return (p1.x-p.x)*(p2.y-p.y)-(p1.y-p.y)*(p2.x-p.x);
}
node zero(0,0);
struct poly {
    int n;
    node p[510];
    void read() {
        for(int i=0; i<n; i++) p[i].read();
    }
    double area() {
        double ans=0;
        for(int i=0; i<n; i++) ans+=mult(zero,p[i],p[(i+1)%n]);
        return fabs(ans/2);
    }
};
bool ins(node p1,node p2,node p3,node p4,node &tmp) {
    double a,b,c,a1,b1,c1;
    a=p1.y-p2.y;
    b=p2.x-p1.x;
    c=p2.y*p1.x-p1.y*p2.x;
    a1=p3.y-p4.y;
    b1=p4.x-p3.x;
    c1=p4.y*p3.x-p3.y*p4.x;
    if(sgn(a*b1-b*a1)==0) return 0;
    tmp=node((b*c1-b1*c)/(a*b1-b*a1),(a*c1-a1*c)/(a1*b-a*b1));
    return 1;
}
poly rec;
void cut(node p1,node p2) {
    int i,j;
    poly temp;
    node tmp;
    double a,b;
    temp.n=0;
    for(i=0; i<rec.n; i++) {
        j=(i+1)%rec.n;
        a=mult(p1,p2,rec.p[i]);
        b=mult(p1,p2,rec.p[j]);
        if(sgn(a)<=0) temp.p[temp.n++]=rec.p[i];
        if(sgn(a*b)<0) {
            ins(p1,p2,rec.p[i],rec.p[j],tmp);
            temp.p[temp.n++]=tmp;
        }
    }
    rec=temp;
}
double run(node p1,node p2,node p3,node p4) {
    int a=sgn(mult(zero,p1,p2)),b=sgn(mult(zero,p3,p4));
    int oper=a*b;
    if(oper==0) return 0;
    rec.n=3;
    rec.p[0]=zero;
    rec.p[1]=p1;
    rec.p[2]=p2;
    if(a>0) swap(rec.p[1],rec.p[2]);
    if(b>0) {
        cut(zero,p4);
        cut(p4,p3);
        cut(p3,zero);
    } else {
        cut(zero,p3);
        cut(p3,p4);
        cut(p4,zero);
    }
    return oper*rec.area();
}
poly A,B;
void solve() {
    int i,j;
    double ans=A.area()+B.area(),s=0;
    for(i=0; i<A.n; i++)
        for(j=0; j<B.n; j++) {
            s+=run(A.p[i],A.p[(i+1)%A.n],B.p[j],B.p[(j+1)%B.n]);
        }
    printf("%.2lf\n",ans-fabs(s));
}
int main() {
    scanf("%d%d",&A.n,&B.n);
    A.read();
    B.read();
    solve();
}
//费马点
#include<iostream>
#include<math.h>
using namespace std;
const double eps=1e-8;
struct node {
    double x,y;
    node() {}
    node(double _x,double _y):x(_x),y(_y) {}
} p[40];
int n;
double get_dis(node pp) {
    double ans=0;
    int i;
    for(i=0; i<n; i++) ans+=sqrt((pp.x-p[i].x)*(pp.x-p[i].x)+(pp.y-p[i].y)*(pp.y-p[i].y));
    return ans;
}
int main() {
    int t,i;
    node st,temp,nt;
    double tmp,tt,step;
    bool ok;
    scanf("%d",&t);
    while(t--) {
        scanf("%d",&n);
        for(i=0; i<n; i++) scanf("%lf%lf",&p[i].x,&p[i].y);
        st=p[0];
        step=2000;
        tmp=get_dis(st);
        while(step>eps) {
            ok=1;
            while(ok) {
                ok=0,nt=st;
                temp=node(st.x,st.y+step);
                tt=get_dis(temp);
                if(tt<tmp) tmp=tt,ok=1,nt=temp;
                temp=node(st.x,st.y-step);
                tt=get_dis(temp);
                if(tt<tmp) tmp=tt,ok=1,nt=temp;
                temp=node(st.x+step,st.y);
                tt=get_dis(temp);
                if(tt<tmp) tmp=tt,ok=1,nt=temp;
                temp=node(st.x-step,st.y);
                tt=get_dis(temp);
                if(tt<tmp) tmp=tt,ok=1,nt=temp;
                st=nt;
            }
            step/=2.0;
        }
        printf("%.6lf\n",tmp);
    }
    return 0;
}
