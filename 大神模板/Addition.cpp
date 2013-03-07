//单纯形   TJU 3723
#include<iostream>
#include<stdio.h>
#include<string.h>
using namespace std;
const double eps = 1e-8;
const double inf = 1e15 ;
#define OPTIMAL -1
#define UNBOUNDED -2
#define FEASIBLE -3
#define INFEASIBLE -4
#define PIVOTOK 1
const int maxn=1000;
const int maxm=1000;
int basic[maxn],row[maxm],col[maxn];
double c0[maxn];
double dcmp(double x) {
    if(x>eps) return 1;
    if(x<-eps) return -1;
    return 0 ;
}
int Pivot(int n ,int m,double * c,double a[maxn][maxn],double* rhs,int &i,int &j) {
    double min = inf ;
    int k = -1;
    for ( j =0; j<=n ; j++) if(!basic[j]&&dcmp(c[j])>0)
            if(k<0||dcmp(c[j]-c[k])>0) k=j ;
    j=k;
    if(k<0)return OPTIMAL;
    for(k=-1,i=1; i<=m; i++) if(dcmp(a[i][j])>0)
            if(dcmp(rhs[i]/a[i][j]-min)<0) {
                min=rhs[i]/a[i][j];
                k=i;
            }
    i=k;
    if(k<0) return UNBOUNDED;
    else return PIVOTOK;
}
int PhaseII(int n,int m,double* c,double a[maxn][maxn],double* rhs,double &ans,int PivotIndex) {
    int i,j,k,l;
    double tmp;
    while(k=Pivot(n,m,c,a,rhs,i,j),k==PIVOTOK||PivotIndex) {
        if(PivotIndex) {
            j=0;
            i=PivotIndex;
            PivotIndex=0;
        }
        basic[row[i]]=0;
        col[row[i]]=0;
        basic[j]=1;
        col[j]=i;
        row[i]=j;
        tmp=a[i][j];
        for(k=0; k<=n; k++)a[i][k]/=tmp;
        rhs[i]/=tmp;
        for(k=1; k<=m; k++)if(k!=i&&dcmp(a[k][j])) {
                tmp=-a[k][j];
                for(l=0; l<=n; l++)a[k][l]+=tmp*a[i][l];
                rhs[k]+=tmp*rhs[i];
            }
        tmp=-c[j];
        for(l=0; l<=n; l++)c[l]+=a[i][l]*tmp;
        ans-=tmp*rhs[i];
    }
    return k ;
}
int PhaseI(int n,int m,double* c,double a[maxn][maxn],double* rhs,double &ans) {
    int i,j,k=-1;
    double tmp,min= 0,ans0= 0;
    for(i=1; i<=m; i++)if(dcmp(rhs[i]-min)<0) {
            min=rhs[i];
            k=i;
        }
    if(k<0) return FEASIBLE;
    for(i=1; i<=m; i++)a[i][0]=-1;
    for(j=1; j<=n; j++)c0[j]=0;
    c0[0]=-1;
    PhaseII(n,m,c0,a,rhs,ans0,k);
    if(dcmp(ans0)<0)return INFEASIBLE;
    for(i=1; i<=m; i++)a[i][0]=0;
    for(j=1; j<=n; j++)if(dcmp(c[j])&&basic[j]) {
            tmp=c[j];
            ans+=rhs[col[j]]*tmp;
            for(i=0; i<=n; i++)c[i]-=tmp*a[col[j]][i];
        }
    return FEASIBLE;
}
int simplex(int n,int m,double * c,double a[maxn][maxn],double* rhs,double &ans,double* x) { //standard form
    int i,j,k;
    for(i=1; i<=m; i++) {
        for(j=n+1; j<=n+m; j++) a[i][j]=0;
        a[i][n+i]=1;
        a[i][0]=0;
        row[i]=n+i;
        col[n+i]=i;
    }
    k=PhaseI(n+m,m,c,a,rhs,ans);
    if(k==INFEASIBLE) return k;
    k=PhaseII(n+m,m,c,a,rhs,ans,0);
    for(j=0; j<=n+m; j++) x[j]=0;
    for(i=1; i<=m; i++)x[row[i]]=rhs[i];
    return k;
}
int n,m;
double c[maxn],ans,a[maxm][maxn],rhs[maxm],x[maxn];
void init() {
    memset(c,0,sizeof(c));
    memset(basic,0,sizeof(basic));
}
double tc[maxn],tmp[maxm][maxn],tt[maxm];
int fmax() {
    int i,j;
    for(i=1; i<=n; i++) {
        for(j=1; j<=3; j++) {
            a[i][j]=tmp[i][j];
        }
        rhs[i]=tt[i];
    }
    double ans=0;
    init();
    for(i=1; i<=3; i++) c[i]=tc[i];
    int judge=simplex (3 , n, c , a , rhs , ans , x );
    if(judge==OPTIMAL) return dcmp(ans);
    return 1;
}
int fmin() {
    int i,j;
    for(i=1; i<=n; i++) {
        for(j=1; j<=3; j++)
            a[i][j]=tmp[i][j];
        rhs[i]=tt[i];
    }
    double ans=0;
    init();
    for(i=1; i<=3; i++) c[i]=-tc[i];
    int judge=simplex (3 , n, c , a , rhs , ans , x );
    if(judge==OPTIMAL) return dcmp(0-ans);
    return -1;
}
int main() {
    freopen("1004.in","r",stdin);
    freopen("1.out","w",stdout);
    int T,i,f1,j,f2;
    char str[20];
    double b;
    scanf("%d",&T);
    while(T--) {
        scanf("%d%d",&n,&m);
        for(i=1; i<=n; i++) {
            scanf("%s",str);
            for(j=1; j<=3; j++) {
                scanf("%lf",&tmp[i][j]);
            }
            for(j=1; j<=3; j++) {
                scanf("%lf",&b);
                if(str[0]=='B') tmp[i][j]=b-tmp[i][j];
                else tmp[i][j]-=b;
            }
            tt[i]=-eps;
        }
        for(i=n+1; i<=n+3; i++) {
            for(j=1; j<=3; j++) {
                tmp[i][j]=0;
                if(i-n==j) tmp[i][j]=1;
            }
            tt[i]=100;
        }
        n+=3;
        for(i=n+1; i<=n+3; i++) {
            for(j=1; j<=3; j++) {
                tmp[i][j]=0;
                if(i-n==j) tmp[i][j]=-1;
            }
            tt[i]=-1;
        }
        n+=3;
        while(m--) {
            for(i=1; i<=3; i++) scanf("%lf",&tc[i]);
            for(i=1; i<=3; i++) {
                scanf("%lf",&b);
                tc[i]-=b;
            }
            f1=fmax();
            f2=fmin();
            if(f1*f2<=0) printf("I don't know\n");
            else if(f1>0) printf("Burning\n");
            else printf("AngleLover\n");
        }
    }
}
//左偏树 zju 2334
#include<iostream>
using namespace std;
const int N=100010;
int fa[N];
inline int root(int a) {
    if(fa[a]==a)return a;
    return fa[a]=root(fa[a]);
}
inline void merge(int a,int b) {
    a=root(a);
    b=root(b);
    fa[a]=b;
}
struct node {
    int k,dis;
    int l,r,f;
    node() {}
    node(int a):k(a),dis(0),l(0),r(0),f(0) {}
} tr[2*N];
inline int troot(int a) {
    if(!a)return 0;
    while(tr[a].f)a=tr[a].f;
    return a;
}
inline int tmerge(int x,int y) {
    if(!(x*y))return x+y;
    if(tr[x].k<tr[y].k)swap(x,y);
    int r=tmerge(tr[x].r,y);
    tr[x].r=r;
    tr[r].f=x;
    if(tr[r].dis>tr[tr[x].l].dis)swap(tr[x].l,tr[x].r);
    if(tr[x].r) tr[x].dis=tr[tr[x].r].dis+1;
    else tr[x].dis=0;
    return x;
}
inline int tins(int i,int k,int& root) {
    tr[i]=node(k);
    return root=tmerge(root,i);
}
inline node tpop(int& root) {
    node out=tr[root];
    int l=tr[root].l,r=tr[root].r;
    tr[root]=node(tr[root].k);
    tr[l].f=tr[r].f=0;
    root=tmerge(l,r);
    return out;
}
inline int tdel(int u) {
    if(!u)return u;
    int x,y,l,r;
    l=tr[u].l;
    r=tr[u].r;
    y=tr[u].f;
    tr[u]=node(tr[u].k);
    tr[x=tmerge(l,r)].f=y;
    if(y&&tr[y].l==u)tr[y].l=x;
    if(y&&tr[y].r==u)tr[y].r=x;
    for(; y; x=y,y=tr[y].f) {
        if(tr[tr[y].r].dis>tr[tr[y].l].dis)swap(tr[y].l,tr[y].r);
        if(tr[tr[y].r].dis+1==tr[y].dis)break;
        tr[y].dis=tr[tr[y].r].dis+1;
    }
    if(x)return troot(x);
    else return troot(y);
}
inline int tadd(int i,int v) { //tr[i].key+=v
    if(!i)return i;
    if(!(tr[i].l+tr[i].r+tr[i].f)) {
        tr[i].k+=v;
        return i;
    }
    int kk=tr[i].k+v;
    int rt=tdel(i);
    return tins(i,kk,rt);
}
int main() {
    int i,n,ask,a,b;
    while(~scanf("%d",&n)) {
        for(i=1; i<=n; i++) {
            scanf("%d",&a);
            tr[i]=node(a);
            fa[i]=i;
        }
        scanf("%d",&ask);
        while(ask--) {
            scanf("%d%d",&a,&b);
            if(root(a)==root(b))printf("-1\n");
            else {
                merge(a,b);
                int aa=troot(a),bb=troot(b);
                aa=tadd(aa,(tr[aa].k/2)-tr[aa].k);
                bb=tadd(bb,(tr[bb].k/2)-tr[bb].k);
                aa=tmerge(aa,bb);
                printf("%d\n",tr[aa].k);
            }
        }
    }
}
/*
计数类问题的解法其实都比较“规矩”，基本就是逐位逼近。首先要求：
f[i][j]表示i位数，各位和是j的数的个数，可以通过简单递推求得。
calc(x,k)表示小于等于x的数里面，各位和是k的数的个数，借助f[][]函数可以逐步逼近的求得。
第二问：累加calc(N,i)直到刚刚大于K，这个i就是第K名的各位数字和，然后枚举每一位的数，利用第一问的函数：对于枚举出来的数字，当作位数和是i，计算排在它前面的数的个数，也是一位一位的算。
详细见09年刘聪的《浅谈数位类统计问题》
*/
#include<iostream>
using namespace std;
typedef unsigned long long ll;
ll f[20][200],ten[25];
ll min(ll a,ll b) {
    return a<b?a:b;
}
void init() {
    int i,j,k;
    memset(f,0,sizeof(f));
    f[0][0]=1;
    for(i=1; i<=19; i++) {
        f[i][0]=1;
        for(j=1; j<200; j++) {
            f[i][j]=f[i-1][j];
            for(k=1; k<=9; k++) {
                if(j-k<0) break;
                f[i][j]+=f[i-1][j-k];
            }
        }
    }
    ten[0]=1;
    for(i=1; i<=20; i++) ten[i]=ten[i-1]*10;
}
ll K,N;
int sum(ll x) {
    int ans=0;
    while(x>0) {
        ans+=x%10;
        x/=10;
    }
    return ans;
}
ll calc(ll x,int k) { //<=x,和是k的个数
    if(x==0) return 0;
    int i,j,cnt=1,vv=sum(x),k2=k;
    ll ans=0;
    for(j=19; j>=0; j--) {
        ll num=x/ten[j];
        x%=ten[j];
        for(i=0; i<num; i++) {
            if(k-i<0) break;
            ans+=f[j][k-i];
        }
        k-=num;
        if(k<0) break;
    }
    if(k2==vv) ++ans;
    return ans;
}
int len(ll x) {
    int ans=0;
    while(x>0) {
        x/=10;
        ans++;
    }
    return ans;
}
ll play1(ll x,int v) {
    int i,j,k,w=0;
    ll temp=x;
    w=len(x);
    ll ans=0;
    int ln=len(N);
    for(i=1; i<v; i++) ans+=calc(N,i);
    if(x==0) return ans;
    temp=x*10;
    for(i=w+1; i<=ln; i++) { //补0
        ans+=calc(min(N,temp-1),v)-f[i-1][v];
        temp*=10;
    }
    temp=x;
    for(i=w; i>=1; i--) { //删尾
        ans+=calc(temp,v)-f[i-1][v];
        temp/=10;
    }
    return ans;
}
ll play2(ll x) {
    int i,j,k,has;
    ll xx=x;
    for(i=1;; i++) {
        ll tt=calc(N,i);
        if(x>tt) x-=tt;
        else break;
    }
    has=i; //位数和
    int ln=len(N);
    ll num=0;
    while(play1(num,sum(num))!=xx) {
        num*=10;
        for(i=0; i<=9; i++)
            if(play1(num+i,has)>=xx) break; //是has,而不仅仅是sum(num+i)
        num+=i;
        if(play1(num,sum(num))==xx) break;
        num--;
    }
    return num;
}
int main() {
    int i,j,k;
    init();
    while(1) {
        cin>>N>>K;
        if(N==0&&K==0) break;
        ll one=play1(K,sum(K));
        ll two=play2(K);
        cout<<one<<" "<<two<<endl;
    }
    return 0;
}
//Treap
#include<iostream>
#include<map>
using namespace std;
const int MAX=100000+50;
struct tnode {
    int l,r;
    int v,p;
    int size,weight;
} tn[MAX];
int neww;//1起始
int newtnode(int v) {
    tn[neww].v=v;
    tn[neww].p=rand();
    tn[neww].l=tn[neww].r=0;
    tn[neww].size=tn[neww++].weight=1;
    return neww-1;
}
struct treap {
    int root;
    treap() {
        root=0;
        srand(16);
    }
    void lr(int& x) {
        int y=tn[x].l;
        tn[x].l=tn[y].r;
        tn[y].r=x;
        tn[x].size=tn[tn[x].l].size+tn[tn[x].r].size+tn[x].weight;
        tn[y].size=tn[tn[y].l].size+tn[tn[y].r].size+tn[y].weight;
        x=y;
    }
    void rr(int& x) {
        int y=tn[x].r;
        tn[x].r=tn[y].l;
        tn[y].l=x;
        tn[x].size=tn[tn[x].l].size+tn[tn[x].r].size+tn[x].weight;
        tn[y].size=tn[tn[y].l].size+tn[tn[y].r].size+tn[y].weight;
        x=y;
    }
    bool contain(int x,int v) {
        if(!x)return false;
        if(tn[x].v==v)return true;
        else if(tn[x].v<v)return contain(tn[x].r,v);
        else return contain(tn[x].l,v);
    }
    void insert(int& x,int v) {
        if(!x) {
            x=newtnode(v);
            return;
        }
        tn[x].size++;
        if(tn[x].v>v) {
            insert(tn[x].l,v);
            if(tn[tn[x].l].p<tn[x].p) lr(x);
        } else if(tn[x].v<v) {
            insert(tn[x].r,v);
            if(tn[tn[x].r].p<tn[x].p)rr(x);
        } else tn[x].weight++;
    }
    void remove(int& x,int v) {
        if(x==root&&(!contain(x,v)))return;
        tn[x].size--;
        if(tn[x].v>v) remove(tn[x].l,v);
        else if(tn[x].v<v) remove(tn[x].r,v);
        else if((--tn[x].weight)==0)
            erase(x);
    }
    void erase(int &x) {
        if(tn[x].l+tn[x].r==0) x=0;
        else if(tn[x].l*tn[x].r==0)x=tn[x].l+tn[x].r;
        else if(tn[tn[x].l].p<tn[tn[x].r].p) {
            lr(x);
            erase(tn[x].r);
        } else {
            rr(x);
            erase(tn[x].l);
        }
    }
    int find_kth(int x,int v) {
        if(v>tn[x].size)return -1;
        if(v<=tn[tn[x].l].size)return find_kth(tn[x].l,v);
        else if (tn[tn[x].l].size+tn[x].weight<v)
            return find_kth(tn[x].r,v-tn[tn[x].l].size-tn[x].weight);
        return tn[x].v;
    }
    int rank(int x,int v) {
        if(!x)return 0;
        if (v<tn[x].v) return rank(tn[x].l,v);
        else return tn[tn[x].l].size+tn[x].weight+rank(tn[x].r,v);
    }
} ta,tb;
map<int,int>M;
int inv[MAX],offset,n;
void top(int v) {
    if (M.count(v)>0)ta.remove(ta.root,M[v]);
    M[v]=--offset;
    ta.insert(ta.root,M[v]);
    inv[MAX-1-offset]=v;
    if (!tb.contain(tb.root,v))tb.insert(tb.root,v);
}
int query(int v) {
    if(M.count(v)>0)return ta.rank(ta.root,M[v]);
    else return v+tn[tb.root].size-tb.rank(tb.root,v);
}
int rank(int v) {
    if(v<=tn[ta.root].size)
        return inv[MAX-1-ta.find_kth(ta.root,v)];
    else if(tn[ta.root].size==0)return v;
    else {
        int l,r,mid;
        v-=tn[ta.root].size;
        for(l=1,r=n; l<=r;) {
            mid=(l+r)>>1;
            if(mid-tb.rank(tb.root,mid)>=v)r=mid-1;
            else l=mid+1;
        }
        return l;
    }
}
int main() {
    int TC,tc,Q,v;
    char ask[30];
    tn[0].l=tn[0].r=tn[0].size=tn[0].weight=0;
    scanf("%d",&TC);
    for(tc=1; tc<=TC; tc++) {
        neww=1;
        offset=MAX;
        ta.root=tb.root=0;
        M.clear();
        scanf("%d%d",&n,&Q);
        printf("Case %d:\n",tc);
        while(Q--) {
            scanf("%s%d",ask,&v);
            switch(ask[0]) {
            case 'T':
                top(v);
                break;
            case 'R':
                printf("%d\n",rank(v));
                break;
            case 'Q':
                printf("%d\n",query(v));
                break;
            };
        }
    }
}
//2进制插头
#include<iostream>
using namespace std;
int map[15][15];
__int64 dp[15][15][1<<12];
int main() {
    int TC,tc=1,r,c,i,j,k,lim;
    scanf("%d",&TC);
    while(TC--) {
        scanf("%d%d",&r,&c);
        for(i=1; i<=r; i++)for(j=1; j<=c; j++)
                scanf("%d",&map[i][j]);

        dp[0][c][0]=1;
        for(i=1,lim=(1<<(c+1)); i<=r; i++) {
            for(j=0; j<(lim>>1); j++)
                dp[i][0][j<<1]=dp[i-1][c][j];
            for(j=1; j<=c; j++)for(k=0; k<lim; k++) {
                    int p=1<<j;
                    int q=p>>1;
                    bool x=k&p;
                    bool y=k&q;
                    if(map[i][j]) {
                        dp[i][j][k]=dp[i][j-1][k^p^q];
                        if(x^y)dp[i][j][k]+=dp[i][j-1][k];
                    } else {
                        if(x|y)dp[i][j][k]=0;
                        else dp[i][j][k]=dp[i][j-1][k];
                    }
                }
        }
        printf("Case %d: There are %I64d ways to eat the trees.\n",tc++,dp[r][c][0]);
    }
}
//3进制插头
#include<iostream>
using namespace std;
const int HA=100009;
const int MN=15;
bool map[MN][MN];//地图，判断线段可放不
int jz[MN],n,m,nn,mm;
__int64 ans;//jz: 进制移动的位数
int myhash[HA];//myhash里存的东东是当前状态在state里第几个位置
int now;//now：滚动用
__int64 sum[2][100000];//state状态，sum:该状态的和
int state[2][100000];
int tot[2];//该格子状态的总数
void Hash_in(int s,__int64 data) {
    int hashpos=s%HA;
    while (myhash[hashpos]) {
        if (state[now][myhash[hashpos]]==s) { //如果状态存在
            sum[now][myhash[hashpos]]+=data;
            return;
        }
        hashpos++;
        if (hashpos==HA) hashpos=0;
    }
    myhash[hashpos]=tot[now];
    state[now][tot[now]]=s;
    sum[now][tot[now]]=data;
    tot[now]++;
}
void solove() {
    int i,j,temps;
    ans=0;
    now=0;
    tot[now]=1;
    state[0][0]=0;
    sum[0][0]=1;

    for(i=1; i<=n; i++) {
        for(j=1; j<=m; j++) {
            now^=1;
            tot[now]=0;
            memset(myhash,0,sizeof(myhash));
            memset(sum[now],0,sizeof(sum[now]));
            memset(state[now],0,sizeof(state[now]));

            for(int k=0; k<tot[now^1]; k++) {
                __int64 s=state[now^1][k],ss=sum[now^1][k];
                int p=(s>>jz[j-1])%4;//取第j-1位的括号
                int q=(s>>jz[j])%4;//取第j位的括号

                if(!map[i][j]) { //map[i][j]不可放,00->00
                    if(!(p+q))Hash_in(s,ss);
                } else if(!(p+q)) { //00->()
                    if(map[i][j+1]&&map[i+1][j]) {
                        temps=s+1*(1<<jz[j-1])+2*(1<<jz[j]);
                        Hash_in(temps,ss);
                    }
                } else if((!p)&&q) { //0(->0(,0(->(0,0)->0),0)->)0
                    if(map[i][j+1])Hash_in(s,ss);
                    if(map[i+1][j]) {
                        temps=s+q*(1<<jz[j-1])-q*(1<<jz[j]);
                        Hash_in(temps,ss);
                    }
                } else if(p&&(!q)) {
                    if(map[i+1][j])Hash_in(s,ss);
                    if(map[i][j+1]) {
                        temps=s-p*(1<<jz[j-1])+p*(1<<jz[j]);
                        Hash_in(temps,ss);
                    }
                } else if(p+q==2) {
                    temps=s-1*(1<<jz[j-1])-1*(1<<jz[j]);
                    int mark=1;//寻找匹配的括号
                    for(int ii=j+1; ii<=m; ii++) {
                        int u=(s>>jz[ii])%4;
                        if(u==1)mark++;
                        else if(u==2)mark--;
                        if(!mark) {
                            temps=temps-2*(1<<jz[ii])+1*(1<<jz[ii]);
                            break;
                        }
                    }
                    Hash_in(temps,ss);
                } else if(p+q==4) {
                    temps=s-2*(1<<jz[j-1])-2*(1<<jz[j]);
                    int mark=1;//寻找匹配的括号
                    for(int ii=j-2; ii>=1; ii--) {
                        int u=(s>>jz[ii])%4;
                        if(u==2)mark++;
                        else if(u==1)mark--;
                        if(!mark) {
                            temps=temps-1*(1<<jz[ii])+2*(1<<jz[ii]);
                            break;
                        }
                    }
                    Hash_in(temps,ss);
                } else if(p==2) {
                    temps=s-2*(1<<jz[j-1])-1*(1<<jz[j]);
                    Hash_in(temps,ss);
                } else if(q==2) {
                    if(i==nn&&j==mm)ans+=ss;
                }
            }
        }
        for (int ii=0; ii<tot[now]; ii++) //每一行做完后，应该乘以4
            state[now][ii]=((state[now][ii]<<2)&((1<<(2*(m+1)))-1));
    }
}
int main() {
    int i,j;
    for(i=0; i<MN; i++)jz[i]=i<<1;
    while(~scanf("%d%d",&n,&m)) {
        for(i=1; i<=n; i++) {
            char tem[20];
            scanf("%s",tem);
            for(j=1; j<=m; j++) {
                map[i][j]=(tem[j-1]=='.');
                if(map[i][j])nn=i,mm=j;
            }
        }

        solove();
        printf("%I64d\n",ans);
    }
}
//4进制独立插头
#include<iostream>
using namespace std;
const int HA=10009;
const int MN=10;
int map[MN][MN],jz[MN],n,m;
int myhash[HA];//myhash里存的东东是当前状态在state里第几个位置
int ans,now,state[2][10000],sum[2][10000];
int tot[2];//该格子状态的总数
void Hash_in(int s,int data) {
    int hashpos=s%HA;
    while (myhash[hashpos]) {
        if (state[now][myhash[hashpos]]==s) { //如果状态存在
            if(sum[now][myhash[hashpos]]<data)sum[now][myhash[hashpos]]=data;
            return;
        }
        hashpos++;
        if (hashpos==HA) hashpos=0;
    }
    myhash[hashpos]=tot[now];
    state[now][tot[now]]=s;
    sum[now][tot[now]]=data;
    tot[now]++;
}
int getk(int st,int j,int p) {
    int tcnt=1,t,k,di,dd;
    if(p==1) {
        di=1 ;
        dd=2;
    } else {
        di=-1;
        dd=1;
    }
    for(k=j;; k+=di) {
        t=(st>>(k<<1))&3;
        if(t==p)tcnt++;
        else if(t==dd) {
            tcnt--;
            if(tcnt==0)return k;
        }
    }
}
void bug(int a) {
    printf(" ");
    while(a) {
        printf("%d",a%2);
        a/=2;
    }
    printf(" ");
}
void solove() {
    int i,j,ii,temps;
    now=0;
    tot[now]=1;
    state[now][0]=sum[now][0]=0;
    for(i=1; i<=n; i++) {
        for(j=1; j<=m; j++) {
            now^=1;
            tot[now]=0;
            memset(myhash,0,sizeof(myhash));
            memset(sum[now],0,sizeof(sum[now]));
            memset(state[now],0,sizeof(state[now]));

            for(int k=0; k<tot[now^1]; k++) {
                int s=state[now^1][k],ss=sum[now^1][k]+map[i][j];
                int p=(s>>jz[j-1])&3;//取第j-1位的括号
                int q=(s>>jz[j])&3;//取第j位的括号
                int c3=(s>>jz[m+1])&3;//独立插头个数
                //bug(s);
                //printf("%d %d %d %d %d\n",i,j,p,q,c3);
                if(!(p+q)) { // 当前格子为空
                    Hash_in(s,ss-map[i][j]);// 可以不走这个格子
                    if(!map[i][j])continue;

                    if(i!=n&&j!=m) { //00->()
                        temps=s+1*(1<<jz[j-1])+2*(1<<jz[j]);
                        Hash_in(temps,ss);
                    }
                    if(c3<2) {
                        if(j!=m) {
                            temps=s+3*(1<<jz[j])+1*(1<<jz[m+1]);//增加向右的独立插头
                            //	printf("!%d\n",temps);
                            Hash_in(temps,ss);
                        }
                        if(i!=n) {
                            temps=s+3*(1<<jz[j-1])+1*(1<<jz[m+1]);// 增加向下的独立插头
                            Hash_in(temps,ss);
                        }
                    }
                } else if(p&&q) {
                    if(!map[i][j])continue;
                    if(p==q) {
                        if(p==3) {
                            if(ans<ss)ans=ss;    // 两个为独立插头，直接连接，更新答案
                        } else if(p==1) {
                            temps=s-1*(1<<jz[j-1])-1*(1<<jz[j]);
                            int ii=getk(s,j+1,1);
                            temps=temps-2*(1<<jz[ii])+1*(1<<jz[ii]);
                            Hash_in(temps,ss);
                        } else {
                            temps=s-2*(1<<jz[j-1])-2*(1<<jz[j]);
                            int ii=getk(s,j-2,2);
                            temps=temps-1*(1<<jz[ii])+2*(1<<jz[ii]);
                            Hash_in(temps,ss);
                        }
                    } else if(p==3||q==3) {
                        temps=s-p*(1<<jz[j-1])-q*(1<<jz[j]);
                        if(p==3)ii=getk(s,q==1?j+1:j-2,q);
                        else ii=getk(s,p==1?j+1:j-2,p);
                        temps|=(3<<jz[ii]);
                        Hash_in(temps,ss);
                    } else if(p==2) {
                        temps=s-2*(1<<jz[j-1])-1*(1<<jz[j]);
                        Hash_in(temps,ss);
                    } else {
                        if(i==n&&j==m)Hash_in(0,ss);
                    }
                } else {
                    if(p) {
                        if(p==3) {
                            //	printf("!\n");
                            if(c3==1&&(s&((1<<jz[m+1])-1-(3<<jz[j-1])))==0&&ans<ss) {
                                // 独立插头为一端，并且其他都是连通，更新答案
                                ans=ss;
                                //	printf("!!\n");
                            }
                        } else if(c3<2) { //把括号另一端变成独立插头，要把这个插头封住。
                            temps=s-(p<<jz[j-1])+(1<<jz[m+1]);
                            ii=getk(s,p==1?j+1:j-2,p);
                            temps|=(3<<jz[ii]);
                            Hash_in(temps,ss);
                        }
                        if(map[i][j]==0) continue; // 注意要在前面转移完之后判断。
                        if(j!=m) {
                            temps=s-(p<<jz[j-1])+(p<<jz[j]);
                            Hash_in(temps,ss);// 直接向右连通
                        }
                        if(i!=n) Hash_in(s,ss); // 直接向下连通
                    }
                    if(q) {
                        if(q==3) {
                            if(c3==1&&(s&((1<<jz[m+1])-1-(3<<jz[j])))==0&&ans<ss)
                                // 独立插头为一端，并且其他都是连通，更新答案
                                ans=ss;
                        } else if(c3<2) { //把括号另一端变成独立插头，要把这个插头封住。
                            temps=s-(q<<jz[j])+(1<<jz[m+1]);
                            ii=getk(s,q==1?j+1:j-2,q);
                            temps|=(3<<jz[ii]);
                            Hash_in(temps,ss);
                        }
                        if(map[i][j]==0) continue; // 注意要在前面转移完之后判断。
                        if(j!=m)Hash_in(s,ss);// 直接向右连通
                        if(i!=n) { // 直接向下连通
                            temps=s+(q<<jz[j-1])-(q<<jz[j]);
                            Hash_in(temps,ss);
                        }
                    }
                }
            }
        }
        for(ii=0; ii<tot[now]; ii++) { //每一行做完后，应该乘以4
            int c3=(state[now][ii]>>jz[m+1])&3;//独立插头个数
            state[now][ii]=((state[now][ii]<<2)&((1<<(2*(m+1)))-1));
            state[now][ii]|=c3<<jz[m+1];
        }
    }
}
int main() {
    int TC,i,j;
    freopen("data.in","r",stdin);
    freopen("data.out","w",stdout);
    for(i=0; i<MN; i++)jz[i]=i<<1;
    scanf("%d",&TC);
    while(TC--) {
        memset(map,0,sizeof(map));
        scanf("%d%d",&n,&m);
        for(ans=0,i=1; i<=n; i++)for(j=1; j<=m; j++) {
                scanf("%d",&map[i][j]);
                if(map[i][j]>ans)ans=map[i][j];
            }
        solove();
        printf("%d\n",ans);
    }
}
/*
pku 1741
1） 每次算出该子树的重心，以该重心进行遍历；
2） 求出以该点为根的所有链路长度，然后排序后o(n）得出 d[a] + d[b] <= k的对数
3） 减去重复的d[a] + d[b] <= k；
*/
#include<iostream>
#include<algorithm>
using namespace std;
#define FOR(a) for(int i=map[a];~i;i=stor[i].next)
const int MN=10001;
struct node {
    int v,w,next;
    node() {}
    node(int a,int b,int c):v(a),w(b),next(c) {}
} stor[2*MN];
int neww,map[MN],sum[MN],bal[MN],sz,now[MN],n,lim,D[MN],ans;
bool use[MN];
void add(int s,int t,int v) {
    stor[neww]=node(t,v,map[s]);
    map[s]=neww++;
    stor[neww]=node(s,v,map[t]);
    map[t]=neww++;
}
void dfs(int u,int fa) {
    sum[u]=bal[u]=0;
    FOR(u)if(stor[i].v!=fa&&(!use[stor[i].v])) {
        dfs(stor[i].v,u);
        sum[u]+=sum[stor[i].v];
        bal[u]=max(bal[u],sum[stor[i].v]);
    }
    sum[u]++;
    now[sz++]=u;
}
int getroot(int u) {
    int tot,min,i;
    sz=0;
    dfs(u,-1);
    tot=sum[u];
    for(min=n,u=-1,i=0; i<sz; i++)
        if(min>max(bal[now[i]],tot-sum[now[i]])) {
            min=max(bal[now[i]],tot-sum[now[i]]);
            u=now[i];
        }
    return u;
}
void getdeepth(int u,int de,int fa) {
    D[sz++]=de;
    FOR(u)if(stor[i].v!=fa&&(!use[stor[i].v]))
        getdeepth(stor[i].v,de+stor[i].w,u);
}
void calc(int t) {
    int l,r;
    sort(D,D+sz);
    for(l=0,r=sz-1; l<r;)
        if(D[l]+D[r]<=lim) {
            ans+=t*(r-l);
            l++;
        } else r--;
}
void solove(int root) {
    root=getroot(root);
    sz=0;
    getdeepth(root,0,-1);
    calc(1);
    use[root]=1;
    FOR(root)if(!use[stor[i].v]) {
        sz=0;
        getdeepth(stor[i].v,stor[i].w,root);
        calc(-1);
    }
    FOR(root)if(!use[stor[i].v])
        solove(stor[i].v);
}
int main() {
    int i;
    while(scanf("%d%d",&n,&lim),n+lim) {
        neww=ans=0;
        for(i=1; i<=n; i++)map[i]=-1,use[i]=0;
        for(i=1; i<n; i++) {
            int s,t,v;
            scanf("%d%d%d",&s,&t,&v);
            add(s,t,v);
        }
        solove(1);
        printf("%d\n",ans);
    }
}
/*
spoj 375
[问题简述]
给定一棵带边权的树,要求设计一个算法,使之能够完成以下2种操作:
1.[Query]给定v1,v2,求出v1到v2路径上的边权的最大值
2.[Change]修改某一条边的边权
LCA使用倍增算法
有了LCA,我们只要求一个节点到他的一个祖先的路径的边权最大值就行了.
**将轻边归入重链处理.
我不知道大家是否认为这是不值得一提的技巧,至少我在膜拜oimaster代码之前不是这样做的(真是太傻×了).
其实,如果把轻边归入重链,那么每条长度为n的重链便由n条边(1条轻边,(n-1)条重边)和n个节点组成.
这样每条边都属于1条重链.
经过根节点的重链是没有轻边的,我们可以加一条虚拟的,边权为-inf的边,以便统一处理.
 */
#include<iostream>
using namespace std;
#define FOR(x) for(int i=map[x];~i;i=stor[i].next)
#define Tree(a) root[belong[a]], 1, 1, length[belong[a]]
#define Left fro,(u<<1),l,mid
#define Right fro,(u<<1)+1,mid+1,r
#define si stor[i]
const int MN=10010;
const int Log=15;
const int inf=0x3fffffff;
struct node {
    int v,w,id,next;
    node() {}
    node(int a,int b,int c,int d):v(a),w(b),id(c),next(d) {}
} stor[3*MN];
int neww,map[MN],weight[MN],depth[MN],sum[MN],n;
int tmax[5*MN],root[5*MN],op,id2v[MN];
int belong[MN],pos[MN],length[MN],anc[MN][Log];
void add(int s,int t,int i,int w) {
    stor[neww]=node(t,w,i,map[s]);
    map[s]=neww++;
    stor[neww]=node(s,w,i,map[t]);
    map[t]=neww++;
}
int NewST(int len) {
    op+=4*len;
    return op-4*len;
}
void fix(int fro,int u,int l,int r,int p,int w) { //fro树根,u偏移量
    if(l==r) {
        tmax[fro+u]=w;
        return;
    }
    int mid=(l+r)>>1;
    if(p>mid)fix(Right,p,w);
    else fix(Left,p,w);
    tmax[fro+u]=max(tmax[fro+(u<<1)],tmax[fro+1+(u<<1)]);
}
int getmax(int fro,int u,int l,int r,int ll,int rr) { //fro树根,u偏移量
    if(ll<=l&&rr>=r)return tmax[fro+u];
    int mid=(l+r)>>1;
    if(rr<=mid)return getmax(Left,ll,rr);
    if(ll>mid)return getmax(Right,ll,rr);
    return max(getmax(Left,ll,mid),getmax(Right,mid+1,rr));
}
void dfs1(int u,int fa) {
    sum[u]=1;
    anc[u][0]=fa;
    depth[u]=depth[fa]+1;
    for(int i=1; i<Log&&anc[u][i-1]; i++)
        anc[u][i]=anc[anc[u][i-1]][i-1];
    FOR(u)if(si.v!=fa) {
        weight[si.v]=si.w;//
        id2v[si.id]=si.v;//
        dfs1(si.v,u);
        sum[u]+=sum[si.v];
    }
}
void dfs2(int u,int len,int fa) {
    int ma=0,tem=0;
    FOR(u)if(si.v!=fa&&sum[si.v]>ma)
        ma=sum[si.v],tem=si.v;
    if(!tem) {
        int p,i;
        for(p=u,i=len-1; i; i--)p=anc[p][0];
        root[p]=NewST(len);
        length[p]=len;
        for(i=len; i; i--,u=anc[u][0]) {
            belong[u]=p;
            pos[u]=i;
            fix(Tree(u),i,weight[u]);
        }
        return;
    }
    dfs2(tem,len+1,u);
    FOR(u)if(si.v!=fa&&si.v!=tem)
        dfs2(si.v,1,u);
}
void BuildTree() {
    memset(anc,0,sizeof(anc));
    depth[0]=0;
    weight[1]=0-inf;
    dfs1(1,0);
    memset(tmax,0,sizeof(tmax));//-inf
    dfs2(1,1,0);
}
int LCA(int a,int b) {
    if(depth[a]>depth[b])swap(a,b);
    for(int t=depth[b]-depth[a],k=0; k<Log; k++)
        if(t&(1<<k))b=anc[b][k];
    if(a==b)return a;
    for(int k=Log-1; k>=0; k--)if(anc[a][k]!=anc[b][k])
            a=anc[a][k],b=anc[b][k];
    return anc[a][0];
}
void Change(int id,int w) {
    weight[id2v[id]]=w;
    fix(Tree(id2v[id]),pos[id2v[id]],w);
}
inline void Up(int &a, int b) {
    if (a < b) a = b;
}
int Query(int a,int b) {
    int res=-inf;
    while(a!=b) {
        if(belong[a]!=belong[b]) {
            Up(res,getmax(Tree(b),1,pos[b]));
            b=anc[belong[b]][0];
        } else {
            Up(res,getmax(Tree(b),pos[a]+1,pos[b]));
            break;
        }
    }
    return res;
}
void solove() {
    int s,t,w,i;
    scanf("%d",&n);
    neww=op=1;
    for(i=1; i<=n; i++)map[i]=-1;
    for(i=1; i<n; i++) {
        scanf("%d%d%d",&s,&t,&w);
        add(s,t,i,w);
    }
    BuildTree();
    static char cmd[100];
    while(scanf("%s",cmd)) {
        if(cmd[0]=='D')break;
        scanf("%d%d",&s,&t);
        if(cmd[0]=='Q') {
            int c=LCA(s,t);
            printf("%d\n",max(Query(c,s),Query(c,t)));
        } else Change(s,t);
    }
}
int main() {
    int TC,flag;
    scanf("%d",&TC);
    flag=0;
    while(TC--) {
        if(!flag)flag=1;
        else printf("\n");
        solove();
    }
}
//spoj  1825
#include<iostream>
#include<algorithm>
using namespace std;
#define FOR(a) for(int i=map[a];~i;i=stor[i].next)
const int inf=0x3fffffff;
const int MN=200001;
struct node {
    int v,w,next;
    node() {}
    node(int a,int b,int c):v(a),w(b),next(c) {}
} stor[2*MN];
struct tt {
    int v,w,val;
    tt() {}
    tt(int a,int b,int c):v(a),w(b),val(c) {}
    bool operator <(const tt& x)const {
        return w<x.w;
    }
} que[MN];
int neww,n,lim,ans,map[MN],sum[MN],bal[MN],bs[MN],now[MN],sz,g[MN],f[MN];
bool cro[MN],use[MN];
void add(int s,int t,int v) {
    stor[neww]=node(t,v,map[s]);
    map[s]=neww++;
    stor[neww]=node(s,v,map[t]);
    map[t]=neww++;
}
void dfs(int u,int fa) {
    sum[u]=bal[u]=bs[u]=0;
    FOR(u)if(stor[i].v!=fa&&(!use[stor[i].v])) {
        dfs(stor[i].v,u);
        sum[u]+=sum[stor[i].v];
        bs[u]+=bs[stor[i].v];
        bal[u]=max(bal[u],sum[stor[i].v]);
    }
    sum[u]++;
    bs[u]+=cro[u];
    now[sz++]=u;
}
int getroot(int u) {
    int tot,min,i;
    sz=0;
    dfs(u,-1);
    tot=sum[u];
    for(min=n,u=-1,i=0; i<sz; i++)
        if(min>max(bal[now[i]],tot-sum[now[i]])) {
            min=max(bal[now[i]],tot-sum[now[i]]);
            u=now[i];
        }
    return u;
}
void getg(int u,int sb,int sw,int fa) {
    g[sb]=max(g[sb],sw);
    FOR(u)if(stor[i].v!=fa&&(!use[stor[i].v]))
        getg(stor[i].v,sb+cro[stor[i].v],sw+stor[i].w,u);
}
void solove(int root,int fa,int fro) {
    int to=fro,j;
    root=getroot(root);
    use[root]=1;
    FOR(root)if(stor[i].v!=fa&&(!use[stor[i].v])) {
        que[to++]=tt(stor[i].v,bs[stor[i].v],stor[i].w);
        solove(stor[i].v,root,to);
    }
    sort(que+fro,que+to);
    int limm=lim-cro[root];
    int ft=-1;
    for(int i=fro; i<to; i++)if(!use[que[i].v]) {
            for(j=0; j<=que[i].w; j++)g[j]=0-inf;
            getg(que[i].v,cro[que[i].v],que[i].val,root);
            int ma=0-inf;
            if(i==fro) {
                for(j=0; j<=que[i].w&&j<=limm; j++)
                    ma=max(ma,g[j]),f[j]=ma;
                ft=min(limm,que[i].w);
            } else {
                for(j=0; j<=que[i].w&&j<=limm; j++) {
                    int tem=limm-j;
                    if(tem>ft)tem=ft;
                    if(f[tem]==0-inf||g[j]==0-inf)continue;
                    ans=max(ans,f[tem]+g[j]);
                }
                for(j=0; j<=que[i].w&&j<=limm; j++) {
                    ma=max(ma,g[j]);
                    if(ft>=j) ma=max(ma,f[j]);
                    f[j]=ma;
                }
                ft=min(limm,que[i].w);
            }
        }
    if(ft>=0) ans=max(ans,f[min(ft,limm)]);
    use[root]=0;
}
int main() {
    int i,s,t,v,m;
    while(~scanf("%d%d%d",&n,&lim,&m)) {
        neww=0;
        ans=0;
        for(i=1; i<=n; i++)map[i]=-1,cro[i]=0,use[i]=0;
        while(m--)scanf("%d",&s),cro[s]=1;
        for(i=1; i<n; i++) {
            scanf("%d%d%d",&s,&t,&v);
            add(s,t,v);
        }
        solove(1,-1,0);
        printf("%d\n",ans);
    }
}
//spoj  2666
#include<iostream>
#include<math.h>
using namespace std;
#define FOR(x) for(edge *i=map[x];i;i=i->next)
const int N = 100010,MaxLevel=40,inf=0x3fffffff;
struct edge {
    int s,v,w;
    bool use;
    edge *next,*ver;
    edge() {}
    edge(int a,int b,int c,edge* e,edge* f):
        s(a),v(b),w(c),use(0),next(e),ver(f) {}
} stor[2*N],*ed,*map[N];
void add(int s,int t,int w) {
    (*ed)=edge(s,t,w,map[s],ed+1);
    map[s]=ed++;
    (*ed)=edge(t,s,w,map[t],ed-1);
    map[t]=ed++;
}
int dist[MaxLevel][N],*td;
bool is_left[MaxLevel][N],in_heap[MaxLevel][N];
bool color[N];
struct par {
    int f,s;
    par() {}
    par(int a,int b):f(a),s(b) {}
};
par tmem[N*MaxLevel],*pmem=tmem;
#define top() a[1]
#define cmp(x,y) (x.f<y.f)
struct priority_queue {
    par *a;
    int n;
    void tmalloc(int len) {
        a=pmem;
        pmem+=len+1;
        n=0;
    }
    void push(par t) {
        int i=++n;
        while (i>1 && cmp(a[i>>1],t)) a[i]=a[i>>1], i>>=1;
        a[i]=t;
    }
    void pop() {
        par t=a[n--];
        int i=1, j=2;
        while (j<=n) {
            if (j<n && cmp(a[j],a[j+1])) ++j;
            if (cmp(t,a[j])) a[i]=a[j], i=j, j<<=1;
            else break;
        }
        a[i]=t;
    }
};
struct BT {
    int level,max;
    edge *e;
    BT *lc,*rc;
    priority_queue l,r;
    void dfs(int u,int fa,int lr) {
        is_left[level][u]=lr;
        in_heap[level][u]=1;
        if(lr)l.push(par(dist[level][u],u));
        else r.push(par(dist[level][u],u));
        FOR(u)if(!i->use&&i->v!=fa) {
            dist[level][i->v]=dist[level][u]+i->w;
            dfs(i->v,u,lr);
        }
    }
    void update() {
        while(l.n&&color[l.top().s])in_heap[level][l.top().s]=0,l.pop();
        while(r.n&&color[r.top().s])in_heap[level][r.top().s]=0,r.pop();
        if(l.n&&r.n)max=l.top().f+r.top().f+e->w;
        else max=-inf;
        if(lc&&lc->max>max) max=lc->max;
        if(rc&&rc->max>max) max=rc->max;
    }
    void init(edge *ee,int sz,int lev,int lsize,int rsize) {
        e=ee;
        level=lev;
        lc=rc=NULL;
        l.tmalloc(lsize);
        r.tmalloc(rsize);
        e->use=e->ver->use=1;
        dist[lev][e->s]=dist[lev][e->v]=0;
        dfs(e->s,e->v,1);
        dfs(e->v,e->s,0);
    }
    void remove(int u) {
        if(is_left[level][u]) {
            if(lc)lc->remove(u);
        } else {
            if(rc)rc->remove(u);
        }
        update();
    }
    void insert(int u) {
        if(!in_heap[level][u]) {
            in_heap[level][u]=1;
            if(is_left[level][u])l.push(par(dist[level][u],u));
            else r.push(par(dist[level][u],u));
        }
        if(is_left[level][u]) {
            if(lc)lc->insert(u);
        } else {
            if(rc)rc->insert(u);
        }
        update();
    }
} bt[5*N],*op;
edge* best;
int size[N],minsize;
void dfs1(int u,int fa,int tot) {
    size[u]=1;
    FOR(u)if((!i->use)&&i->v!=fa) {
        dfs1(i->v,u,tot);
        int temp=size[i->v]*2-tot;
        temp=abs(temp);
        if (temp<minsize) minsize=temp, best=i;
        size[u]+=size[i->v];
    }
}
BT* Build(int u,int sz,int lev) {
    BT* now=op++;
    minsize=sz;
    dfs1(u,0,sz);
    int rsize=size[best->v],lsize=sz-size[best->v];
    now->init(best,sz,lev,lsize,rsize);
    if(now->l.n>1)
        now->lc=Build(now->e->s,lsize,lev+1);
    if(now->r.n>1)
        now->rc=Build(now->e->v,rsize,lev+1);
    now->update();
    return now;
}
int main() {
    char ch[3];
    int i,t,ask,n,s,w,white;
    while(~scanf("%d",&n)) {
        for(i=1; i<=n; i++)map[i]=NULL;
        ed=stor;
        op=bt;
        pmem=tmem;
        for(i=1; i<n; i++) {
            scanf("%d%d%d",&s,&t,&w);
            add(s,t,w);
        }
        memset(color,0,sizeof(color));
        white=n;//0为白
        BT* root=Build(1,n,0);

        scanf("%d",&ask);
        while (ask--) {
            scanf("%s",ch);
            if (ch[0]=='A') {
                if (!white) printf("They have disappeared.\n");
                else printf("%d\n",max(0,root->max));
            } else {
                scanf("%d",&t);
                color[t]=!color[t];
                if (color[t]) --white,root->remove(t);
                else ++white,root->insert(t);
            }
        }
    }
}
//hdu 3601
#include<iostream>
using namespace std;
#define FOR(x) for(int i=map[x];~i;i=stor[i].next)
#define Tree(a) root[belong[a]], 1, 1, length[belong[a]]
#define Left fro,(u<<1),l,mid
#define Right fro,(u<<1)+1,mid+1,r
#define si stor[i]
const int MN=30010;
const int Log=17;
const int inf=0x3fffffff;
struct node {
    int v,next;
    node() {}
    node(int a,int b):v(a),next(b) {}
} stor[3*MN];
int neww,map[MN],weight[MN],depth[MN],sum[MN],n;
int tmax[5*MN],tsum[5*MN],root[5*MN],op;
int belong[MN],pos[MN],length[MN],anc[MN][Log];
void add(int a,int b) {
    stor[neww]=node(b,map[a]);
    map[a]=neww++;
    stor[neww]=node(a,map[b]);
    map[b]=neww++;
}
void dfs1(int u,int fa) {
    sum[u]=1;
    anc[u][0]=fa;
    depth[u]=depth[fa]+1;
    for(int i=1; i<Log&&anc[u][i-1]; i++)
        anc[u][i]=anc[anc[u][i-1]][i-1];
    FOR(u)if(si.v!=fa) {
        dfs1(si.v,u);
        sum[u]+=sum[si.v];
    }
}
int NewST(int len) {
    op+=4*len;
    return op-4*len;
}
void dfs2(int u,int len,int fa) {
    int ma=0,tem=0;
    FOR(u)if(si.v!=fa&&sum[si.v]>ma)
        ma=sum[si.v],tem=si.v;
    if(!tem) {
        int p,i;
        for(p=u,i=len-1; i; i--)p=anc[p][0];
        root[p]=NewST(len);
        length[p]=len;
        for(i=len; i; i--,u=anc[u][0]) {
            belong[u]=p;
            pos[u]=i;
        }
        return;
    }
    dfs2(tem,len+1,u);
    FOR(u)if(si.v!=fa&&si.v!=tem)
        dfs2(si.v,1,u);
}
void build() {
    memset(anc,0,sizeof(anc));
    depth[0]=0;
    dfs1(1,0);
    memset(tmax,0,sizeof(tmax));
    memset(tsum,0,sizeof(tsum));
    dfs2(1,1,0);
}
int LCA(int a,int b) {
    if(depth[a]>depth[b])swap(a,b);
    for(int t=depth[b]-depth[a],k=0; k<Log; k++)
        if(t&(1<<k))b=anc[b][k];
    if(a==b)return a;
    for(int k=Log-1; k>=0; k--)if(anc[a][k]!=anc[b][k])
            a=anc[a][k],b=anc[b][k];
    return anc[a][0];
}
void fix1(int fro,int u,int l,int r,int p,int w) { //fro树根,u偏移量
    if(l==r) {
        tmax[fro+u]=w;
        return;
    }
    int mid=(l+r)>>1;
    if(p>mid)fix1(Right,p,w);
    else fix1(Left,p,w);
    tmax[fro+u]=max(tmax[fro+(u<<1)],tmax[fro+1+(u<<1)]);
}
void fix2(int fro,int u,int l,int r,int p,int w) { //fro树根,u偏移量
    if(l==r) {
        tsum[fro+u]+=w;
        return;
    }
    int mid=(l+r)>>1;
    if(p>mid)fix2(Right,p,w);
    else fix2(Left,p,w);
    tsum[fro+u]=tsum[fro+(u<<1)]+tsum[fro+1+(u<<1)];
}
void change(int id,int w) {
    weight[id]+=w;
    fix1(Tree(id),pos[id],weight[id]);
    fix2(Tree(id),pos[id],w);
}
int getmax(int fro,int u,int l,int r,int ll,int rr) { //fro树根,u偏移量
    if(ll<=l&&rr>=r)return tmax[fro+u];
    int mid=(l+r)>>1;
    if(rr<=mid)return getmax(Left,ll,rr);
    if(ll>mid)return getmax(Right,ll,rr);
    return max(getmax(Left,ll,mid),getmax(Right,mid+1,rr));
}
inline void Up(int &a, int b) {
    if (a < b) a = b;
}
int qmax(int a,int b) {
    int res=weight[a];
    while(a!=b) {
        if(belong[a]!=belong[b]) {
            Up(res,getmax(Tree(b),1,pos[b]));
            b=anc[belong[b]][0];
        } else {
            Up(res,getmax(Tree(b),pos[a],pos[b]));
            break;
        }
    }
    return res;
}
int getsum(int fro,int u,int l,int r,int ll,int rr) { //fro树根,u偏移量
    if(ll<=l&&rr>=r)return tsum[fro+u];
    int mid=(l+r)>>1;
    if(rr<=mid)return getsum(Left,ll,rr);
    if(ll>mid)return getsum(Right,ll,rr);
    return getsum(Left,ll,mid)+getsum(Right,mid+1,rr);
}
int qsum(int a,int b) {
    int res=weight[a];
    while(a!=b) {
        if(belong[a]!=belong[b]) {
            res+=getsum(Tree(b),1,pos[b]);
            b=anc[belong[b]][0];
        } else {
            res+=getsum(Tree(b),pos[a]+1,pos[b]);
            break;
        }
    }
    return res;
}
void Print(int fro,int u,int l,int r) {
    printf("%d ",tsum[fro+u]);
    if(l==r)return;
    int mid=(l+r)>>1;
    Print(Right);
    Print(Left);
}
int main() {
    int i,a,b,m,c;
    char ask[30];
    while(~scanf("%d",&n)) {
        neww=op=1;
        for(i=1; i<=n; i++)map[i]=-1,weight[i]=0;
        for(i=1; i<n; i++) {
            scanf("%d%d",&a,&b);
            add(a,b);
        }
        build();
        for(i=1; i<=n; i++) {
            scanf("%d",&a);
            change(i,a);
        }
        scanf("%d",&m);
        while(m--) {
            scanf("%s%d%d",ask,&a,&b);
            if(ask[0]=='C')
                change(a,b-weight[a]);
            else if(ask[1]=='M') {
                c=LCA(a,b);
                printf("%d\n",max(qmax(c,a),qmax(c,b)));
            } else {
                c=LCA(a,b);
                printf("%d\n",qsum(c,a)+qsum(c,b)-weight[c]);
            }
        }
    }
}
//最大团 Hdu 1530
#include<stdio.h>
#define MAX_SIZE 58
int dp[MAX_SIZE];//dp[]
bool inset[MAX_SIZE];
bool map[MAX_SIZE][MAX_SIZE];//map[][]
int N,dfs;
bool findIt;
void Memcpy(bool*d,bool*s) { //memcpy
    for(int i=0; i<N; i++)
        d[i]=s[i];
}
int FindAdj(int start,bool*s) { //find neighbor
    for(int i=start+1; i<N; i++)
        if(s[i])
            return i;
    return -1;
}
void DFS(bool*s,int start,int depth) { //DFS
    int first=FindAdj(start,s),i;
    if(first==-1) {
        if(depth>dfs) {
            dfs=depth;
            findIt=true;
        }
        return;
    }
    bool news[MAX_SIZE];
    while(first!=-1) {
        if(depth+N-start<=dfs||depth+dp[first]<=dfs)
            return;//cut if can not find a best ans
        news[first]=false;
        for(i=first+1; i<N; i++)
            news[i]=s[i]&map[first][i];
        DFS(news,first,depth+1);
        if(findIt)return;
        first=FindAdj(first,s);
    }
}
int main() {
    int i,j;
    while(scanf("%d",&N)&&N) {
        for(dfs=i=0; i<N; i++)
            for(j=0; j<N; j++)
                scanf("%d",&map[i][j]);
        for(i=N-1; i>=0; i--) {
            Memcpy(inset,map[i]);
            findIt=false;
            DFS(inset,i,1);
            dp[i]=dfs;
        }
        printf("%d\n",dp[0]);
    }
    return 0;
}
//HK
#include <stdio.h>
int head,tail,total,n,edgeN;
int q[600],slice[600],lead[600],headg[600],tailg[600];
struct stct0 {
    int next,v;
} edge[260000];
int DFS(int point) {
    int ptrg;
    ptrg=headg[point];
    while(ptrg+1) {
        if(lead[edge[ptrg].v]==-1||(slice[lead[edge[ptrg].v]]==slice[point]+1&&DFS(lead[edge[ptrg].v]))) {
            lead[edge[ptrg].v]=point;
            return 1;
        }
        ptrg=edge[ptrg].next;
    }
    slice[point]=-1;
    return 0;
}
int BFS() {
    int flag,ptrg,i;
    head=0;
    for(i=0; i<=n-1; i++)
        slice[i]=1;
    for(i=0; i<=n-1; i++)
        if(lead[i]!=-1)
            slice[lead[i]]=-1;
    tail=0;
    for(i=0; i<=n-1; i++)
        if(slice[i]==1)
            q[tail++]=i;
    flag=0;
    while(head!=tail) {
        ptrg=headg[q[head]];
        while(ptrg!=-1) {
            if(lead[edge[ptrg].v]==-1)
                flag=1;
            else if(slice[lead[edge[ptrg].v]]==-1) {
                slice[lead[edge[ptrg].v]]=slice[q[head]]+1;
                q[tail++]=lead[edge[ptrg].v];
            }
            ptrg=edge[ptrg].next;
        }
        ++head;
    }
    return flag;
}
void addedge(int a,int b) {
    if(headg[a]==-1)
        headg[a]=edgeN;
    else
        edge[tailg[a]].next=edgeN;
    edge[edgeN].next=-1;
    edge[edgeN].v=b;
    tailg[a]=edgeN++;
    return;
}
main() {
    int p,t,i;
    while(~scanf("%d",&n)) {
        n/=2;
        //init{
        edgeN=0;
        for(i=0; i<=n-1; i++) {
            headg[i]=-1;
            tailg[i]=-1;
        }
        //}
        for(i=0; i<=n-1; i++) {
            scanf("%d",&t);
            while(t--) {
                scanf("%d",&p);
                addedge(i,p-1);
            }
        }
        total=0;
        for(i=0; i<=n-1; i++)
            lead[i]=-1;
        while(BFS()) {
            i=0;
            while(i!=tail&&slice[q[i]]!=2) {
                if(DFS(q[i]))
                    total++;
                i++;
            }
        }
        printf("%d\n",total*2);
    }
    return 0;
}
/*
//带花树
Ural-1099 Work scheduling 难度：*
http://acm.timus.ru/problem.aspx?space=1&num=1099
大意：有N个士兵，每个士兵可以另一个士兵巡逻，给出M组士兵愿意在一起巡逻，求出最大的士兵配置方案。
分析：直接求最大匹配就好了，然后输出所有的组合
ZJU 3316 Game 难度：**
http://acm.zju.edu.cn/onlinejudge/showProblem.do?problemCode=3316
大意：有N个棋子在棋盘上，2个人轮流拿走一个棋子，第一步可以拿任意一个，而之后每一步必须拿上一步拿走的棋子曼哈顿长度L以内的棋子，问，后手是否能 赢
分析：把每一个棋子与周围距离为L的棋子都连上边后，形成一些联通块。易知，一轮游戏只能在某一个联通块里开始直到结束。那么，如果有一个联通块不是完美 匹配，先手就可以走那个没被匹配到的点，后手不论怎么走，都必然走到一个被匹配的点上，先手就可以顺着这个交错路走下去，最后一定是后手没有路可走，因为 如果还有路可走，这一条交错路，就是一个增广路，必然有更大的匹配。
总结：判断是否是完美匹配。
HDU 3446 daizhenyang’s chess 难度：***
http://acm.hdu.edu.cn/showproblem.php?pid=3446
大意：有一个R*C的棋盘，棋盘上有一些格子是幸运的格子，棋盘上只有一个棋子：king。king有一定的走路范围：20个格子如题中图所示。两个人轮 流走，每个人可以移动king走到没被走过的幸运的格子上。问，先手是否能赢。
分析：这题与上一题某方面是一样的。都是与增广路有关。把所有能连的边先连上，我们先不算king，求一次匹配。然后我们再算上king，求一次匹配。如 果第二次匹配的结果比第一次大，说至少明存在一条增广路，且增广路的起点为king所在的点，那么，我们沿着这条路走下去，最后后手必然无路可走。
总结：2次匹配。
HDU 3551 hard problem 难度：****
http://acm.hdu.edu.cn/showproblem.php?pid=3551
大意：给出一个无向图，并给出一个子图的每个点的度数，问是否能去掉一些边，得到这个子图。
分析：因为我们要求出这样的一个子图，是删去一些边，使其度数减少，所以，我们先把每条边拆成2个点，(u,e)和(e,v)，这两个点连一条边，若匹配 这条边表示边e(u,v)不被删除，然后我们把每个点拆成deg[i]-D[i]个点（deg[i]为原图度数，D[i]为子图度数），然后对所有该点和 他的边(i,e)组成的点连一条边，表示这个点缺少的度数，可以由这儿得到。接着，求匹配，如果是完美匹配，即：每个度数都找到了可以使他减少的边。这个 子图可以得到。
总结：让每个需要减少的度数匹配到该点所连接的一条边。
下面是我自己用的带花树的模板，暂时没有出现问题：
补充：
对于2和3的自我理解：
1.假如AB下棋，A先下，那么如果在联通块内有完美匹配，则如果A下了，B必然可以下（仅仅需要走到匹配的点就可以了）！因此如果有完全匹配，B一定可以赢（对应A一定输）。否则反之。
2.对于3题：没有必要像作者说的那样，求两次匹配，只需求一次不带King的匹配就可以了！原理如 [补充1]
*/
struct Graph {
    int n, match[maxn];
    bool adj[maxn][maxn];
    void clear() {
        memset(adj, 0, sizeof(adj));
        n = 0;
    }
    void insert(const int &u, const int &v) {
        get_max(n, max(u, v) + 1);
        adj[u][v] = adj[v][u] = 1;
    }
    int max_match() {
        memset(match, -1, sizeof(match));
        int ans = 0;
        for (int i = 0; i < n; ++i) {
            if (match[i] == -1) {
                ans += bfs(i);
            }
        }
        return ans;
    }
    int Q[maxn], pre[maxn], base[maxn];
    bool hash[maxn];
    bool in_blossom[maxn];
    int bfs(int p) {
        memset(pre, -1, sizeof(pre));
        memset(hash, 0, sizeof(hash));
        for (int i = 0; i < n; ++i) {
            base[i] = i;
        }
        Q[0] = p;
        hash[p] = 1;
        for (int s = 0, t = 1; s < t; ++s) {
            int u = Q[s];
            for (int v = 0; v < n; ++v) {
                if (adj[u][v] && base[u] != base[v] && v != match[u]) {
                    if (v == p || (match[v] != -1 && pre[match[v]] != -1)) {
                        int b = contract(u, v);
                        for (int i = 0; i < n; ++i) {
                            if (in_blossom[base[i]]) {
                                base[i] = b;
                                if (hash[i] == 0) {
                                    hash[i] = 1;
                                    Q[t++] = i;
                                }
                            }
                        }
                    } else if (pre[v] == -1) {
                        pre[v] = u;
                        if (match[v] == -1) {
                            argument(v);
                            return 1;
                        } else {
                            Q[t++] = match[v];
                            hash[match[v]] = 1;
                        }
                    }
                }
            }
        }
        return 0;
    }
    void argument(int u) {
        while (u != -1) {
            int v = pre[u];
            int k = match[v];
            match[u] = v;
            match[v] = u;
            u = k;
        }
    }
    void change_blossom(int b, int u) {
        while (base[u] != b) {
            int v = match[u];
            in_blossom[base[v]] = in_blossom[base[u]] = true;
            u = pre[v];
            if (base[u] != b) {
                pre[u] = v;
            }
        }
    }
    int contract(int u, int v) {
        memset(in_blossom, 0, sizeof(in_blossom));
        int b = find_base(base[u], base[v]);
        change_blossom(b, u);
        change_blossom(b, v);
        if (base[u] != b) {
            pre[u] = v;
        }
        if (base[v] != b) {
            pre[v] = u;
        }
        return b;
    }
    int find_base(int u, int v) {
        bool in_path[maxn] = {};
        while (true) {
            in_path[u] = true;
            if (match[u] == -1) {
                break;
            }
            u = base[pre[match[u]]];
        }
        while (!in_path[v]) {
            v = base[pre[match[v]]];
        }
        return v;
    }
};
//ISAP模板
const int inf=0x7FFFFFFF;
const int v_size=210;
const int e_size=500;
struct node {
    int p,next,f;
    node() {}
    node(int a,int b,int c):p(a),f(b),next(c) {}
} map[v_size],edge[e_size];
int n,m,op,s,t,ans,lim;
int level[v_size],num[v_size],que[v_size],pre[v_size],start[v_size];
void add(int a,int b,int c) {
    edge[op]=node(b,c,map[a].next);
    map[a].next=op++;
}
void bfs() {
    int head,tail,temp,mid,i,j,tt;
    for(i=s; i<=t; i++)
        level[i]=lim+1;
    memset(num,0,sizeof(num));
    head=0;
    tail=1;
    que[0]=t;
    level[t]=1;
    while(head<tail) {
        temp=tail;
        for(i=head; i<temp; i++) {
            mid=que[i];
            num[level[mid]]++;
            for(j=map[mid].next; j+1; j=edge[j].next) {
                if(!(j&1)) continue;
                tt=edge[j].p;
                if(level[tt]>=lim+1&&edge[j^1].f>0) {
                    level[tt]=level[mid]+1;
                    que[tail++]=tt;
                }
            }
        }
        head=temp;
    }
}
int min(int a,int b) {
    return a<b?a:b;
}
void flow() {
    int temp,fin;
    fin=inf;
    temp=t;
    while(temp!=s) {
        fin=min(fin,edge[pre[temp]].f);
        temp=edge[pre[temp]^1].p;
    }
    ans+=fin;
    temp=t;
    while(temp!=s) {
        edge[pre[temp]].f-=fin;
        edge[pre[temp]^1].f+=fin;
        temp=edge[pre[temp]^1].p;
    }
}
bool repeat(int &cur) {
    int fin=lim,mid,temp;
    mid=map[cur].next;
    while(mid!=-1) {
        if(edge[mid].f>0) fin=min(fin,level[edge[mid].p]);
        mid=edge[mid].next;
    }
    temp=level[cur];
    level[cur]=fin+1;
    num[temp]--;
    num[level[cur]]++;
    if(cur!=s) cur=edge[pre[cur]^1].p;
    return num[temp]>0;
}
void solve() {
    bfs();
    int cur,i,mid,tt;
    ans=0;
    cur=s;
    for(i=s; i<=t; i++)
        start[i]=map[i].next;
    while(level[s]<=lim) {
        mid=start[cur];
        while(mid!=-1) {
            tt=edge[mid].p;
            if(level[tt]+1==level[cur]&&edge[mid].f>0) break;
            mid=edge[mid].next;
        }
        if(mid!=-1) {
            start[cur]=mid;
            pre[edge[mid].p]=mid;
            cur=edge[mid].p;
            if(cur==t) {
                flow();
                cur=s;
            }
        } else {
            start[cur]=map[cur].next;
            if(!repeat(cur)) break;
        }
    }
    printf("%d\n",ans);
}
int main() {
    int i,a,b,c;
    while(~scanf("%d%d",&m,&n)) {
        s=1;
        t=n;
        lim=n;
        for(i=s; i<=t; i++)
            map[i].next=-1;
        op=0;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            add(a,b,c);
            add(b,a,0);
        }
        solve();
    }
}
/* 后缀数组＋栈扫描
POJ3415
 把两个串拼接起来，求k前缀子串数量a，再分别求两个串k前缀子串的数量b, c，然后用a – b – c。
 其中求完拼接串的后缀数组，两个串的后缀数组可以由品拼接串的后缀数组直接拆分得到，不需要重
 新O(N * log N)求。

先帖下求单个串的k前缀子串数量核心算法：

假设有m个子串在heigt数组中连续且值等于b，那么他们所产生的k前缀数量为 C(m, 2) * (b – k + 1)

用一个栈扫描height(简称hgt)数组

规则1、如果hgt大于栈顶元素，入栈，i++，循环continue

规则2、如果等于栈顶元素，不入栈，i++，循环continue

规则3、如果小于栈顶元素，计算i与栈顶元素之间的长度，按照入栈规则，这一段的hgt值都是一样的，
      可以容易地按上一段所说的规则进行统计。这里还要分两种情况

      1、hgt小于栈顶下面的第一个元素，我是先把栈顶元素的hgt值“消”成和其下面的一个元素值
        一样，等于将他们合并，统计这部分的差值（下称相对贡献量），；

      2、hgt大于栈顶下面的第一个元素，把栈顶元素的hgt值“消”和hgt一样，其他同1

      不继续扫描i，循环continue

规则4、hgt小于k且栈中有元素，执行统计，且执行出栈操作，不继续扫描i，循环continue

规则5、hgt小于k且栈中无元素，i++，循环continue
 */
#include<iostream>
#include<stdio.h>
#include<string.h>
using namespace std;
const int MAX=200010;
int rank[MAX],trank[MAX],to[MAX],sa[MAX],tsa[MAX],cnt[MAX],h[MAX],a[MAX],b[MAX];
int n;
char str[MAX],temp[MAX];
bool equal(int tmp[],int a,int b,int k) {
    if(a+k>=n&&b+k>=n) return tmp[a]==tmp[b];
    if(a+k<n&&b+k<n) return tmp[a]==tmp[b]&&tmp[a+k]==tmp[b+k];
    return 0;
}
void creat_suffix_array() {
    int i,k;
    for(i=0; i<=400; i++) cnt[i]=0;
    for(i=0; i<n; i++) cnt[str[i]]++;
    for(i=1; i<=400; i++) cnt[i]+=cnt[i-1];
    for(i=n-1; i>=0; i--) sa[--cnt[str[i]]]=i;
    for(rank[sa[0]]=0,i=1; i<n; i++) {
        rank[sa[i]]=rank[sa[i-1]];
        if(str[sa[i]]!=str[sa[i-1]]) rank[sa[i]]++;
    }
    for(k=1; k<n; k<<=1) {
        for(i=0; i<n; i++) cnt[rank[sa[i]]]=i+1;
        for(i=n-1; i>=0; i--)
            if(sa[i]>=k) tsa[--cnt[rank[sa[i]-k]]]=sa[i]-k;
        for(i=n-k; i<n; i++) tsa[--cnt[rank[i]]]=i;
        for(trank[tsa[0]]=0,i=1; i<n; i++) {
            trank[tsa[i]]=trank[tsa[i-1]];
            if(!equal(rank,tsa[i],tsa[i-1],k)) trank[tsa[i]]++;
        }
        for(i=0; i<n; i++) sa[i]=tsa[i],rank[i]=trank[i];
    }
}
void creat_height() {
    int i,j,k;
    for(k=0,i=0; i<n; i++) {
        if(rank[i]==n-1) h[rank[i]]=k=0;
        else {
            if(k) k--;
            j=sa[rank[i]+1];
            for(; i+k<n&&j+k<n&&str[i+k]==str[j+k]; k++);
            h[rank[i]]=k;
        }
    }
}
int get(int a,int b) {
    int c=a-b+1;
    if(c<0) return 0;
    return c;
}
int main() {
    __int64 ans,ss;
    int k,len,top,i,t;
    while(scanf("%d",&k),k) {
        scanf("%s",str);
        len=strlen(str);
        str[len]='#';
        scanf("%s",str+len+1);
        n=strlen(str);
        creat_suffix_array();
        creat_height();
        for(i=0; i<n; i++) to[i]=sa[i]<len,h[i]=get(h[i],k);
        a[0]=-1,b[0]=0;
        ans=0;
        for(t=0; t<=1; t++)
            for(i=0,ss=0,top=0; i<n-1; i++) {
                a[++top]=h[i];
                b[top]=to[i]==t;
                ss+=(__int64)a[top]*(__int64)b[top];
                while(a[top]<=a[top-1]) {
                    ss-=(__int64)(a[top-1]-a[top])*(__int64)b[top-1];
                    b[top-1]+=b[top];
                    a[top-1]=a[top];
                    top--;
                }
                if(to[i+1]!=t) ans+=ss;
            }
        printf("%I64d\n",ans);
    }
    return 0;
}
//平面最近点对
#include<cstdio>
#include<cstdlib>
#include<cmath>
#define size 1000

using namespace std;

struct point {
    int x,y;
};

point rec[size],xx[size],yy[size];

int cmp(const void *a,const void *b) {
    return ((point *)a)->x-((point *)b)->x;
}

int dis(point p,point q) {
    return (p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y);
}

int solve(point *temp1,point *temp2,int start,int end,int y_pos[]) {
    int mid,min,i,j,t,u,v,t1,t2,l;
    int yr[size],yl[size],k[size];
    point u1,v1,u2,v2;
    bool judge[size];
    if(end-start<=2) {
        min=999999;
        for(i=start; i<=end-1; i++)
            for(j=i+1; j<=end; j++) {
                t=dis(rec[xx[i].y],rec[xx[j].y]);
                if(min>t) {
                    min=t;
                    *temp1=rec[xx[i].y];
                    *temp2=rec[xx[j].y];
                }
            }
        return min;
    }
    mid=(start+end)>>1;
    for(i=start; i<=end; i++)
        judge[xx[i].y]=(i>=mid+1);
    u=mid+1;
    v=start;
    for(i=start; i<=end; i++)
        if(judge[y_pos[i]]) yr[u++]=y_pos[i];
        else yl[v++]=y_pos[i];
    t1=solve(&u1,&v1,start,mid,yl);
    t2=solve(&u2,&v2,mid+1,end,yr);
    if(t1<t2) {
        l=-1;
        mid=(rec[xx[mid].y].x+rec[xx[mid+1].y].x)>>1;
        for(i=start; i<=end; i++)
            if(abs(rec[y_pos[i]].x-mid)<t1) k[++l]=y_pos[i];
        *temp1=u1;
        *temp2=v1;
        for(i=0; i<=l-1; i++)
            for(j=i+1; j<=l&&j<=i+7; j++) {
                if(t1>dis(rec[k[i]],rec[k[j]])) {
                    t1=dis(rec[k[i]],rec[k[j]]);
                    *temp1=rec[k[i]];
                    *temp2=rec[k[j]];
                }
            }
        return t1;
    } else {
        l=-1;
        mid=(rec[xx[mid].y].x+rec[xx[mid+1].y].x)>>1;
        for(i=start; i<=end; i++)
            if(abs(rec[y_pos[i]].x-mid)<t2) k[++l]=y_pos[i];
        *temp1=u2;
        *temp2=v2;
        for(i=0; i<=l-1; i++)
            for(j=i+1; j<=l&&j<=i+7; j++) {
                if(t2>dis(rec[k[i]],rec[k[j]])) {
                    t2=dis(rec[k[i]],rec[k[j]]);
                    *temp1=rec[k[i]];
                    *temp2=rec[k[j]];
                }
            }
        return t2;
    }
}

int main() {
    int n,i;
    int y_pos[size];
    point ans1,ans2;
    scanf("%d",&n);
    for(i=0; i<=n-1; i++) {
        scanf("%d%d",&rec[i].x,&rec[i].y);
        xx[i].x=rec[i].x;
        xx[i].y=i;
        yy[i].x=rec[i].y;
        yy[i].y=i;
    }
    qsort(xx,n,sizeof(point),cmp);
    qsort(yy,n,sizeof(point),cmp);
    for(i=0; i<=n-1; i++)
        y_pos[i]=yy[i].y;
    printf("%.3lf\n",sqrt(solve(&ans1,&ans2,0,n-1,y_pos)));
    printf("(%d,%d),(%d,%d)\n",ans1.x,ans1.y,ans2.x,ans2.y);

    return 0;
}

