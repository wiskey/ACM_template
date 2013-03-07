无向图求割顶
struct node {
    int p,next;
    node(int a=0,int b=0):p(a),next(b) {}
} edge[100010],map[110];
int num[110],deep[110],low[110];
bool ans[110];
int op;
void dfs(int pos,int pre) {
    int i,tt;
    low[pos]=deep[pos]=op++;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        tt=edge[i].p;
        if(tt==pre) continue;
        if(deep[tt]) low[pos]=min(low[pos],deep[tt]);
        else {
            dfs(tt,pos);
            low[pos]=min(low[pos],low[tt]);
            num[pos]++;
            if(pos==1&&num[1]>1) ans[pos]=1;
            if(pos!=1&&low[tt]>=deep[pos]) ans[pos]=1;
        }
    }
}
void add(int a,int b) {
    edge[op]=node(b,map[a].next);
    map[a].next=op++;
}
int main() {
    int n,a,b,sum,i;
    char ch;
    while(scanf("%d",&n),n) {
        for(i=1; i<=n; i++)
            map[i].next=-1;
        op=0;
        while(scanf("%d",&a),a) {
            scanf("%c",&ch);
            if(ch=='\n') continue;
            while(1) {
                scanf("%d%c",&b,&ch);
                add(a,b);
                add(b,a);
                if(ch=='\n') break;
            }
        }
        sum=0;
        op=1;
        memset(num,0,sizeof(num));
        memset(ans,0,sizeof(ans));
        memset(deep,0,sizeof(deep));
        dfs(1,0);
        for(i=1; i<=n; i++)
            if(ans[i]) sum++;
        printf("%d\n",sum);
    }
}
//无向图求桥（边的双连通分量）
struct node {
    int p,next,cover;
    node(int a=0,int b=0,int c=0):p(a),cover(b),next(c) {}
} edge[10010],map[5010];
int n,m,op;
int low[5010],deep[5010],in[5010];
bool mark[10010];
class set {
private:
    int father[5010];
public:
    void reset(int num) {
        int i;
        for(i=1; i<=num; i++)
            father[i]=i;
    }
    int find(int a) {
        if(a==father[a]) return a;
        return father[a]=find(father[a]);
    }
    void merge(int a,int b) {
        int ta=find(a);
        int tb=find(b);
        father[ta]=tb;
    }
    bool check(int a,int b) {
        return find(a)==find(b);
    }
} S;
void add(int a,int b) {
    edge[op]=node(b,1,map[a].next);
    map[a].next=op++;
}
void dfs(int pos,int pre) {
    low[pos]=deep[pos]=op++;
    int i,tt;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        tt=edge[i].p;
        if(tt==pre) continue;
        if(low[tt]) low[pos]=min(low[pos],deep[tt]);
        else {
            dfs(tt,pos);
            low[pos]=min(low[pos],low[tt]);
            if(low[tt]>deep[pos]) mark[i]=mark[i^1]=1;//重边特别处理
        }
    }
}
void solve() {
    memset(mark,0,sizeof(mark));
    memset(in,0,sizeof(in));
    memset(low,0,sizeof(low));
    memset(deep,0,sizeof(deep));
    op=0;
    dfs(1,0);
    S.reset(n);
    int i,a,b,ans;
    for(i=0; i<2*m; i+=2)
        if(mark[i]==0) S.merge(edge[i].p,edge[i^1].p);
    for(i=0; i<2*m; i+=2) {
        a=edge[i].p;
        b=edge[i^1].p;
        if(S.check(a,b)) continue;
        in[S.find(a)]++;
        in[S.find(b)]++;
    }
    ans=0;
    for(i=1; i<=n; i++)
        ans+=(in[i]==1);
    printf("%d\n",(ans+1)>>1);
}
int main() {
    int mm,i,a,b;
    while(~scanf("%d%d",&n,&m)) {
        for(i=1; i<=n; i++)
            map[i].next=-1;
        op=0;
        mm=m;
        while(mm--) {
            scanf("%d%d",&a,&b);
            add(a,b);
            add(b,a);
        }
        solve();
    }
}

//第K短路A＊
const int inf=0x7FFFFFFF;
int n,m,op,_op,s,t,k;
int dis[1010],in[1010],que[2000000];
bool mark[1010];
struct gg {
    int p,len;
    gg() {}
    gg(int a,int b):p(a),len(b) {}
    bool operator<(const gg a)const {
        return dis[p]+len>dis[a.p]+a.len;
    }
};
struct node {
    int p,w,next;
    node() {}
    node(int a,int b,int c):p(a),w(b),next(c) {}
} edge[100010],_edge[100010],map[1010],_map[1010];
void add(int a,int b,int c) {
    edge[op]=node(b,c,map[a].next);
    map[a].next=op++;
}
void _add(int a,int b,int c) {
    _edge[_op]=node(b,c,_map[a].next);
    _map[a].next=_op++;
}
void spfa() {
    int i;
    for(i=1; i<=n; i++) {
        mark[i]=0;
        dis[i]=inf;
    }
    dis[t]=0;
    int j,head(0),tail(1),temp,tt;
    que[0]=t;
    mark[t]=1;
    while(head<tail) {
        temp=que[head++];
        for(j=_map[temp].next; j+1; j=_edge[j].next) {
            tt=_edge[j].p;
            if(dis[tt]>dis[temp]+_edge[j].w) {
                dis[tt]=dis[temp]+_edge[j].w;
                if(!mark[tt]) {
                    mark[tt]=1;
                    que[tail++]=tt;
                }
            }
        }
        mark[temp]=0;
    }
}
void astar() {
    int i,temp,tt,j,ll;
    priority_queue <gg> H;
    H.push(gg(s,0));
    for(i=1; i<=n; i++)
        in[i]=0;
    while(!H.empty()) {
        temp=H.top().p;
        ll=H.top().len;
        H.pop();
        in[temp]++;
        if(in[t]==k) {
            printf("%d\n",ll);
            return ;
        }
        if(in[temp]>k) continue;
        for(j=map[temp].next; j+1; j=edge[j].next) {
            tt=edge[j].p;
            H.push(gg(tt,ll+edge[j].w));
        }
    }
    printf("-1\n");
}
int main() {
    int i,a,b,c;
    while(~scanf("%d%d",&n,&m)) {
        op=0;
        _op=0;
        for(i=1; i<=n; i++)
            map[i].next=_map[i].next=-1;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            add(a,b,c);
            _add(b,a,c);
        }
        scanf("%d%d%d",&s,&t,&k);
        if(s==t) k++;
        spfa();
        astar();
    }
}
次小生成树
struct node {
    int s,t,l;
    bool operator<(const node a)const {
        return l<a.l;
    }
} rec[10010];
struct gg {
    int p,w;
    gg() {}
    gg(int a,int b):p(a),w(b) {}
};
vector <gg> H[110];
int father[110];
bool mark[10010];
int dis[110][110];
int ans,n,m;
int find(int x) {
    if(father[x]==x) return x;
    return father[x]=find(father[x]);
}
bool check(int x,int y) {
    int tx=find(x),ty=find(y);
    return tx==ty;
}
void merge(int x,int y) {
    int tx=find(x),ty=find(y);
    father[tx]=ty;
}
void kruscal() {
    int times(0),i,a,b;
    for(i=1; i<=n; i++) {
        H[i].clear();
        father[i]=i;
    }
    ans=0;
    for(i=0; i<m; i++)
        mark[i]=0;
    sort(rec,rec+m);
    for(i=0; i<m; i++) {
        a=rec[i].s;
        b=rec[i].t;
        if(check(a,b)) continue;
        merge(a,b);
        ans+=rec[i].l;
        mark[i]=1;
        H[a].push_back(gg(b,rec[i].l));
        H[b].push_back(gg(a,rec[i].l));
        times++;
        if(times==n-1) break;
    }
}
void dfs(int start,int pos,int pre,int len) {
    dis[start][pos]=len;
    int tt;
    vector <gg> ::iterator it;
    for(it=H[pos].begin(); it!=H[pos].end(); it++) {
        tt=(*it).p;
        if(tt==pre) continue;
        dfs(start,tt,pos,max(len,(*it).w));
    }
}
void solve() {
    int i,a,b;
    for(i=0; i<m; i++) {
        if(mark[i]) continue;
        a=rec[i].s;
        b=rec[i].t;
        if(rec[i].l==dis[a][b]) {
            printf("Not Unique!\n");
            return ;
        }
    }
    printf("%d\n",ans);
}
int main() {
    int t,i;
    scanf("%d",&t);
    while(t--) {
        scanf("%d%d",&n,&m);
        for(i=0; i<m; i++)
            scanf("%d%d%d",&rec[i].s,&rec[i].t,&rec[i].l);
        kruscal();
        for(i=1; i<=n; i++)
            dfs(i,i,0,0);
        solve();
    }
}
//强连通分量Kosaraju
struct node {
    int p,next;
    node(int a=0,int b=0):p(a),next(b) {}
} edge[2010],map[310];
bool mark[310];
int op,len;
int rec[310];
void add(int a,int b) {
    edge[op]=node(b,map[a].next);
    map[a].next=op++;
}
void dfs_h(int pos) {
    int i,tt;
    mark[pos]=1;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        if(i&1) continue;
        tt=edge[i].p;
        if(mark[tt]) continue;
        dfs_h(tt);
    }
    rec[len++]=pos;
}
void dfs_t(int pos) {
    int i,tt;
    mark[pos]=1;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        if(!(i&1)) continue;
        tt=edge[i].p;
        if(mark[tt]) continue;
        dfs_t(tt);
    }
}
int main() {
    int n,m,i,a,b,ans;
    while(scanf("%d%d",&n,&m),n+m) {
        for(i=1; i<=n; i++)
            map[i].next=-1;
        op=0;
        while(m--) {
            scanf("%d%d",&a,&b);
            add(a,b);
            add(b,a);
        }
        len=0;
        memset(mark,0,sizeof(mark));
        for(i=1; i<=n; i++)
            if(!mark[i]) dfs_h(i);
        memset(mark,0,sizeof(mark));
        ans=0;
        for(i=n; i>=1; i--)
            if(!mark[rec[i]]) {
                ans++;
                dfs_t(rec[i]);
            }
        printf("%d\n",ans);
    }
}

2-SAT问题
基本逻辑关系
x   1 (~x->x)
0 (x->~x)
x&y 1 (~x->x)^(~y->y)
0 (x->~y)^(y->~x)
x|y 1 (~x->y)^(~y->x)
0 (x->~x)^(y->~y)
x^y 1 (x->~y)^(y->~x)^(~x->y)^(~y->x)
0 (x->y)^(~y->~x)^(~x->~y)^(y->x)
struct node {
    int p,next;
    node() {}
    node(int a,int b):p(a),next(b) {}
} edge[10000000],map[2010],gra[2010];
int part[2010],rec[2010],vis[2010],f[2010];
bool mark[2010];
int n,op,len,ans,nn;
void add(int a,int b) {
    edge[op]=node(b,map[a].next);
    map[a].next=op++;
}
void Plus(int a,int b) {
    edge[op]=node(b,gra[a].next);
    gra[a].next=op++;
}
void dfshead(int pos) {
    int i,tt;
    mark[pos]=1;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        if(i&1) continue;
        tt=edge[i].p;
        if(mark[tt]) continue;
        dfshead(tt);
    }
    rec[++len]=pos;
}
void dfstail(int pos) {
    int i,tt;
    mark[pos]=1;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        if(!(i&1)) continue;
        tt=edge[i].p;
        if(mark[tt]) continue;
        dfstail(tt);
    }
    part[pos]=ans;
}
void dfs(int pos) {
    int i,tt;
    mark[pos]=1;
    for(i=gra[pos].next; i+1; i=edge[i].next) {
        tt=edge[i].p;
        if(mark[tt]) continue;
        dfs(tt);
    }
    rec[++len]=pos;
}
void solve() {
    int i,j,tt;
    memset(mark,0,sizeof(mark));
    len=0;
    for(i=0; i<nn; i++)
        if(!mark[i]) dfshead(i);
    memset(mark,0,sizeof(mark));
    ans=0;
    for(i=nn; i>=1; i--) {
        if(!mark[rec[i]]) {
            ans++;
            f[ans]=rec[i];
            dfstail(rec[i]);
        }
    }
    for(i=0; i<n; i++) {
        if(part[i]==part[i+n]) {
            printf("bad luck\n");
            return ;
        }
    }
    for(i=1; i<=ans; i++)
        gra[i].next=-1;
    for(i=0; i<nn; i++) {
        for(j=map[i].next; j+1; j=edge[j].next) {
            if(j&1) continue;
            tt=edge[j].p;
            if(part[i]==part[tt]) continue;
            Plus(part[i],part[tt]);
        }
    }
    memset(mark,0,sizeof(mark));
    len=0;
    for(i=1; i<=ans; i++)
        if(!mark[i]) dfs(i);
    memset(vis,0,sizeof(vis));
    for(i=1; i<=ans; i++) {
        if(vis[rec[i]]) continue;
        vis[rec[i]]=1;
        tt=f[rec[i]];
        tt-=n;
        if(tt<0) tt+=(n<<1);
        vis[part[tt]]=2;
    }
    for(i=1; i<n; i++) {
        if(i!=1) printf(" ");
        printf("%d%c",i,1==vis[part[i]]?'w':'h');
    }
    printf("\n");
}
int main() {
    int to[200];
    int m,i,a,b;
    char cha,chb;
    to['h']=0;
    to['w']=1;
    while(scanf("%d%d",&n,&m),n+m) {
        nn=n<<1;
        for(i=0; i<nn; i++)
            map[i].next=-1;
        op=0;
        while(m--) {
            scanf("%d%c",&a,&cha);
            scanf("%d%c",&b,&chb);
            add(a+n*to[cha],b+n*(to[chb]^1));
            add(b+n*(to[chb]^1),a+n*to[cha]);
            add(b+n*to[chb],a+n*(to[cha]^1));
            add(a+n*(to[cha]^1),b+n*to[chb]);
        }
        add(n,0);
        add(0,n);
        solve();
    }
}
//点的双连通分量＋交叉染色求奇环
struct node {
    int p,next;
    node(int a=0,int b=0):p(a),next(b) {}
} edge[200010],map[100010];
bool ans[100010];
int op,len,n;
int low[100010],deep[100010],color[100010],stack[100010];
bool dfs(int pos,int pre,int col) {
    deep[pos]=low[pos]=op++;
    color[pos]=col;
    stack[++len]=pos;
    bool res=0,judge;
    int i,tt;
    for(i=map[pos].next; i+1; i=edge[i].next) {
        tt=edge[i].p;
        if(tt==pre) continue;
        if(low[tt]) {
            low[pos]=min(low[pos],deep[tt]);
            if(deep[tt]<deep[pos]&&color[tt]==col) res|=1;
        } else {
            judge=dfs(tt,pos,3-col);
            low[pos]=min(low[pos],low[tt]);
            if(low[tt]<deep[pos]) res|=judge;
            if(low[tt]>=deep[pos]) {
                while(stack[len]!=tt) {
                    ans[stack[len]]|=judge;
                    len--;
                }
                ans[stack[len--]]|=judge;
                ans[pos]|=judge;
            }
        }
    }
    return res;
}
void solve() {
    op=1;
    int i,sum=0;
    for(i=1; i<=n; i++) {
        low[i]=deep[i]=0;
        ans[i]=0;
    }
    for(i=1; i<=n; i++)
        if(!low[i]) {
            len=0;
            dfs(i,0,1);
        }
    for(i=1; i<=n; i++)
        if(ans[i]) sum++;
    printf("%d\n",sum);
}
void add(int a,int b) {
    edge[op]=node(b,map[a].next);
    map[a].next=op++;
}
int main() {
    int t,a,b,i,m;
    scanf("%d",&t);
    while(t--) {
        scanf("%d%d",&n,&m);
        op=0;
        for(i=1; i<=n; i++)
            map[i].next=-1;
        while(m--) {
            scanf("%d%d",&a,&b);
            add(a,b);
            add(b,a);
        }
        solve();
    }
}
//度限度最小生成树

const int inf=0x7FFFFFFF;
map <string,int> S;
int lin[30][30];
struct node {
    int s,t,l;
    node() {}
    node(int a,int b,int c):s(a),t(b),l(c) {}
    bool operator<(const node a)const {
        return l<a.l;
    }
} rec[50000],head[30000],dp[30];
int n,esize,hsize,ans,park,k;
int father[30];
int cost[30][2],dis[30][30];
bool vis[30];
int find(int x) {
    if(father[x]==x) return x;
    return father[x]=find(father[x]);
}
bool check(int x,int y) {
    return find(x)==find(y);
}
void merge(int x,int y) {
    father[find(x)]=find(y);
}
void kruscal() {
    int i,j,a,b;
    sort(rec,rec+esize);
    for(i=1; i<=n; i++)
        father[i]=i;
    for(i=1; i<=n; i++)
        for(j=1; j<=n; j++)
            dis[i][j]=inf;
    ans=0;
    for(i=0; i<esize; i++) {
        a=rec[i].s;
        b=rec[i].t;
        if(check(a,b)) continue;
        ans+=rec[i].l;
        merge(a,b);
        dis[a][b]=dis[b][a]=rec[i].l;
    }
}
void dfs(int pos,int pre,int len,int start,int end) {
    int i;
    dp[pos]=node(start,end,len);
    if(pre==park) dp[pos]=node(start,end,0);
    for(i=1; i<=n; i++) {
        if(i==pre) continue;
        if(dis[pos][i]==inf) continue;
        if(pos==park) {
            dfs(i,pos,0,pos,i);
            continue;
        }
        if(len>dis[pos][i]) dfs(i,pos,len,start,end);
        else dfs(i,pos,dis[pos][i],pos,i);
    }
}
void solve() {
    int i,temp,min,tt,j;
    for(i=1; i<=n; i++) {
        cost[i][0]=inf;
        cost[i][1]=-1;
    }
    for(i=0; i<hsize; i++) {
        temp=head[i].s;
        if(temp==park) temp=head[i].t;
        if(cost[find(temp)][0]>head[i].l) {
            cost[find(temp)][0]=head[i].l;
            cost[find(temp)][1]=temp;
        }
    }
    memset(vis,0,sizeof(vis));
    for(i=1; i<=n; i++) {
        if(cost[i][1]==-1) continue;
        ans+=cost[i][0];
        k--;
        vis[cost[i][1]]=1;
        dis[park][cost[i][1]]=cost[i][0];
        dis[cost[i][1]][park]=cost[i][0];
    }
    while(k) {
        dfs(park,0,0,park,park);
        min=inf;
        tt=-1;
        for(i=0; i<hsize; i++) {
            temp=head[i].s;
            if(temp==park) temp=head[i].t;
            if(vis[temp]) continue;
            if(min>head[i].l-dp[temp].l) {
                tt=temp;
                min=head[i].l-dp[temp].l;
            }
        }
        if(tt==-1) break;
        vis[tt]=1;
        if(min>=0) break;
        ans+=min;
        dis[dp[tt].s][dp[tt].t]=inf;
        dis[dp[tt].t][dp[tt].s]=inf;
        dis[park][tt]=dis[tt][park]=min+dp[tt].l;
        k--;
    }
    printf("Total miles driven: %d\n",ans);
}
int main() {
    int m,l,a,b;
    char sa[20],sb[20];
    string ta,tb;
    while(~scanf("%d",&m)) {
        S.clear();
        n=esize=hsize=0;
        while(m--) {
            scanf("%s%s%d",sa,sb,&l);
            ta=sa;
            tb=sb;
            if(S.find(ta)==S.end()) S[ta]=++n;
            if(S.find(tb)==S.end()) S[tb]=++n;
            a=S[ta];
            b=S[tb];
            if(ta=="Park"||tb=="Park") head[hsize++]=node(a,b,l);
            else rec[esize++]=node(a,b,l);
        }
        scanf("%d",&k);
        park=S["Park"];
        kruscal();
        solve();
    }
}
最优比率生成树

二分法
const double eps=1e-4;
const double inf=1e10;
struct node {
    double x,y,h;
} rec[1010];
double cost[1010],lin[1010][1010],graph[1010][1010],w[1010]

[1010];
bool mark[1010];
double ans;
int n;
double dis(node a,node b) {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
void prim(double temp) {
    int i,j,pos;
    double min;
    for(i=1; i<n; i++)
        for(j=i+1; j<=n; j++)
            graph[i][j]=graph[j][i]=w[i][j]-temp*lin[i][j];
    for(i=1; i<=n; i++) {
        mark[i]=0;
        cost[i]=graph[1][i];
    }
    ans=0;
    mark[1]=1;
    for(i=1; i<n; i++) {
        min=inf;
        for(j=1; j<=n; j++)
            if(!mark[j]&&min>cost[j]) {
                min=cost[j];
                pos=j;
            }
        mark[pos]=1;
        ans+=cost[pos];
        for(j=1; j<=n; j++)
            if(!mark[j]&&cost[j]>graph[pos][j])
                cost[j]=graph[pos][j];
    }
}
void deal() {
    int i,j;
    for(i=1; i<n; i++)
        for(j=i+1; j<=n; j++) {
            lin[j][i]=lin[i][j]=dis(rec[i],rec[j]);
            w[j][i]=w[i][j]=fabs(rec[i].h-rec[j].h);
        }
}
int main() {
    double low,high,mid;
    int i;
    while(scanf("%d",&n),n) {
        for(i=1; i<=n; i++)
            scanf("%lf%lf%lf",&rec[i].x,&rec[i].y,&rec[i].h);
        deal();
        low=0.0;
        high=33.0;
        while(1) {
            mid=(low+high)/2;
            prim(mid);
            if(fabs(ans)<eps) break;
            if(ans>eps) low=mid;
            else high=mid;
        }
        printf("%.3lf\n",mid);
    }
}
迭代逼近法
const double eps=1e-4;
const double inf=1e10;
int n;
bool mark[1010];
int pre[1010];
double w[1010][1010],lin[1010][1010],graph[1010][1010],cost

[1010];
struct node {
    double x,y,h;
} rec[1010];
double sumu,sumd;
void prim(double temp) {
    int i,j,pos;
    double min;
    for(i=1; i<n; i++)
        for(j=i+1; j<=n; j++)
            graph[i][j]=graph[j][i]=w[i][j]-temp*lin[i][j];
    for(i=1; i<=n; i++) {
        mark[i]=0;
        cost[i]=graph[1][i];
        pre[i]=1;
    }
    sumu=0;
    sumd=0;
    mark[1]=1;
    for(i=1; i<n; i++) {
        min=inf;
        for(j=1; j<=n; j++)
            if(!mark[j]&&min>cost[j]) {
                min=cost[j];
                pos=j;
            }
        mark[pos]=1;
        sumu+=w[pre[pos]][pos];
        sumd+=lin[pre[pos]][pos];
        for(j=1; j<=n; j++)
            if(!mark[j]&&cost[j]>graph[pos][j]) {
                cost[j]=graph[pos][j];
                pre[j]=pos;
            }
    }
}
double dis(node a,node b) {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
void deal() {
    int i,j;
    for(i=1; i<n; i++)
        for(j=i+1; j<=n; j++) {
            lin[j][i]=lin[i][j]=dis(rec[i],rec[j]);
            w[j][i]=w[i][j]=fabs(rec[i].h-rec[j].h);
        }
}
int main() {
    int i;
    double cur,pre;
    while(scanf("%d",&n),n) {
        for(i=1; i<=n; i++)
            scanf("%lf%lf%lf",&rec[i].x,&rec[i].y,&rec[i].h);
        deal();
        cur=0;
        pre=0;
        while(1) {
            prim(cur);
            cur=sumu/sumd;
            if(fabs(cur-pre)<eps) break;
            pre=cur;
        }
        printf("%.3lf\n",cur);
    }
}
无向图最小割
int lin[510][510];
int v[510],dis[510];
bool mark[510];
int n,ans;
void solve() {
    int i,max,k,pre,j;
    for(i=1; i<n; i++)
        v[i]=i;
    ans=0x7FFFFFFF;
    while(n>1) {
        max=-1;
        for(i=1; i<n; i++) {
            mark[v[i]]=0;
            dis[i]=lin[0][v[i]];
            if(max==-1||dis[max]<dis[i]) max=i;
        }
        for(i=1; i<=n-1; i++) {
            mark[v[max]]=1;
            if(i==n-1) {
                ans=min(ans,dis[max]);
                for(k=0; k<n; k++)
                    lin[v[k]][v[pre]]=(lin[v[pre]][v[k]]+=lin[v

                                                          [k]][v[max]]);
                v[max]=v[n-1];
            }
            pre=max;
            max=-1;
            for(j=1; j<n; j++) {
                if(!mark[v[j]]) {
                    dis[j]+=lin[v[pre]][v[j]];
                    if(max==-1||dis[j]>dis[max])
                        max=j;
                }
            }
        }
        n--;
    }
}
int main() {
    int m,a,b,c,i,j;
    while(~scanf("%d%d",&n,&m)) {
        for(i=0; i<n; i++)
            for(j=0; j<n; j++)
                lin[i][j]=0;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            lin[a][b]+=c;
            lin[b][a]+=c;
        }
        solve();
        printf("%d\n",ans);
    }
}
LCA问题
const int SIZE=10010;
struct node {
    　　　int p,w;
    　　　node() {}
    　　　node(int a,int b):p(a),w(b) {}
};
int res[1000010];
class set {
private:
    int pre[SIZE];
public:
    void reset(int temp) {
        int i;
        for(i=1; i<=temp; i++)
            pre[i]=i;
    }
    int find(int temp) {
        if(pre[temp]==temp) return temp;
        return pre[temp]=find(pre[temp]);
    }
    void merge(int a,int b) {
        int ta=find(a);
        int tb=find(b);
        pre[ta]=tb;
    }
    bool check(int a,int b) {
        return find(a)==find(b);
    }
} S;
vector <node> G[SIZE];
vector <node> Q[SIZE];
int father[SIZE];
bool mark[SIZE];
int dis[SIZE];
int sear(int pos) {
    if(father[pos]==pos) return pos;
    return father[pos]=sear(father[pos]);
}
void uni(int a,int b) {
    int ta=sear(a);
    int tb=sear(b);
    father[tb]=ta;
}
void lca(int pos,int anc,int len) {
    mark[pos]=1;
    S.merge(pos,anc);
    dis[pos]=len;
    int i;
    for(i=0; i<G[pos].size(); i++) {
        if(!mark[G[pos][i].p]) {
            lca(G[pos][i].p,anc,len+G[pos][i].w);
            uni(pos,G[pos][i].p);
        }
    }
    for(i=0; i<Q[pos].size(); i++) {
        if(mark[Q[pos][i].p]) {
            if(!S.check(anc,Q[pos][i].p)) res[Q[pos][i].w]=-1;
            else res[Q[pos][i].w]=dis[pos]+dis[Q[pos][i].p]-2*dis[sear(Q[pos][i].p)];
        }
    }
}
int main() {
    int n,m,c,i,a,b,q;
    while(~scanf("%d%d%d",&n,&m,&q)) {
        for(i=1; i<=n; i++) {
            G[i].clear();
            Q[i].clear();
        }
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            G[a].push_back(node(b,c));
            G[b].push_back(node(a,c));
        }
        for(i=1; i<=q; i++) {
            scanf("%d%d",&a,&b);
            Q[a].push_back(node(b,i));
            Q[b].push_back(node(a,i));
        }
        memset(mark,0,sizeof(mark));
        S.reset(n);
        for(i=1; i<=n; i++)
            father[i]=i;
        for(i=1; i<=n; i++)
            if(!mark[i]) lca(i,i,0);
        for(i=1; i<=q; i++) {
            if(res[i]==-1) printf("Not connected\n");
            else printf("%d\n",res[i]);
        }
    }
}
多重匹配匈牙利变形
int graph[1010][510];
int map[510][1010];
int num[1010],cur[510];
bool mark[510];
int mid,n;
bool find(int pos) {
    int i,j,temp;
    for(i=0; i<num[pos]; i++) {
        temp=graph[pos][i];
        if(mark[temp]) continue;
        mark[temp]=1;
        if(cur[temp]<mid) {
            map[temp][cur[temp]++]=pos;
            return 1;
        }
        for(j=0; j<mid; j++) {
            if(find(map[temp][j])) {
                map[temp][j]=pos;
                return 1;
            }
        }
    }
    return 0;
}
bool check() {
    int i;
    memset(cur,0,sizeof(cur));
    for(i=1; i<=n; i++) {
        memset(mark,0,sizeof(mark));
        if(!find(i)) return 0;
    }
    return 1;
}
void solve() {
    int high,low;
    high=n;
    low=1;
    while(low<=high) {
        mid=(high+low)>>1;
        if(check()) high=mid-1;
        else low=mid+1;
    }
    if(!check()) mid++;
    printf("%d\n",mid);
}
int main() {
    int m,i,a,j;
    char ch;
    char str[20];
    while(scanf("%d%d",&n,&m),n+m) {
        memset(num,0,sizeof(num));
        for(i=1; i<=n; i++) {
            scanf("%s",str);
            while(1) {
                scanf("%d%c",&a,&ch);
                graph[i][num[i]++]=a;
                if(ch=='\n') break;
            }
        }
        solve();
    }
}
树的同构
const int MUL=1991;
const int MOD=9991;
const int MAX=1010;
int root[2];
int n,state;
int que[MAX],hash[MAX],in[MAX];
bool mark[MAX];
vector <int> map[MAX];
bool cmp(const int a,const int b) {
    return hash[a]<hash[b];
}
void topsort() {
    if(n==2) {
        root[0]=1;
        root[1]=2;
        state=1;
        return ;
    }
    int head,tail,i,temp,mid,tt,j;
    head=tail=0;
    memset(mark,0,sizeof(mark));
    for(i=1; i<=n; i++) {
        if(in[i]==1) {
            que[tail++]=i;
            mark[i]=1;
        }
    }
    while(head<tail) {
        if(tail==n-2) state=1;
        if(tail==n-1) state=0;
        temp=tail;
        for(i=head; i<temp; i++) {
            mid=que[i];
            for(j=0; j<map[mid].size(); j++) {
                tt=map[mid][j];
                if(!mark[tt]) {
                    in[tt]--;
                    if(in[tt]==1) {
                        mark[tt]=1;
                        que[tail++]=tt;
                    }
                }
            }
        }
        head=temp;
    }
    if(state==1) {
        root[0]=que[n-1];
        root[1]=que[n-2];
    }
    if(state==0) root[0]=que[n-1];
}
int dfs(int pos,int pre) {
    int i;
    bool judge=0;
    for(i=0; i<map[pos].size(); i++) {
        if(map[pos][i]!=pre) {
            judge=1;
            hash[map[pos][i]]=dfs(map[pos][i],pos);
        }
    }
    if(!judge) return 1;
    sort(map[pos].begin(),map[pos].end(),cmp);
    int val=1908;
    for(i=0; i<map[pos].size(); i++)
        if(map[pos][i]!=pre) val=((val*MUL)^hash[map[pos][i]])%MOD;
    return val;
}
int main() {
    int t,m,i,a,b,num1,num2,j;
    bool judge;
    int hash1[2],hash2[2];
    scanf("%d",&t);
    while(t--) {
        scanf("%d",&n);
        m=n-1;
        for(i=1; i<=n; i++)
            map[i].clear();
        memset(in,0,sizeof(in));
        while(m--) {
            scanf("%d%d",&a,&b);
            in[a]++;
            in[b]++;
            map[a].push_back(b);
            map[b].push_back(a);
        }
        topsort();
        num1=state;
        if(state==0) {
            hash1[0]=dfs(root[0],0);
        } else {
            hash1[0]=dfs(root[0],0);
            hash1[1]=dfs(root[1],0);
        }
        m=n-1;
        for(i=1; i<=n; i++)
            map[i].clear();
        memset(in,0,sizeof(in));
        while(m--) {
            scanf("%d%d",&a,&b);
            in[a]++;
            in[b]++;
            map[a].push_back(b);
            map[b].push_back(a);
        }
        topsort();
        num2=state;
        if(state==0) {
            hash2[0]=dfs(root[0],0);
        } else {
            hash2[0]=dfs(root[0],0);
            hash2[1]=dfs(root[1],0);
        }
        judge=0;
        if(num1==0&&num2==0) {
            if(hash1[0]==hash2[0]) judge=1;
        }
        if(num1==0&&num2==1) {
            for(i=0; i<=1; i++)
                if(hash1[0]==hash2[i]) judge=1;
        }
        if(num1==1&&num2==0) {
            for(i=0; i<=1; i++)
                if(hash1[i]==hash2[0]) judge=1;
        }
        if(num1==1&&num2==1) {
            for(i=0; i<=1; i++)
                for(j=0; j<=1; j++)
                    if(hash1[i]==hash2[j]) judge=1;
        }
        if(judge) printf("same\n");
        else printf("different\n");
    }
}
网络流
EK模板
const int SIZE=210;
const int INF=0x7FFFFFFF;
int op,ans;
int pre[SIZE],que[SIZE];
bool mark[SIZE];
int n;
struct node {
    int p,next;
    int f;
} edge[SIZE*2],map[SIZE];
void add(int u,int v,int f) {
    edge[op].p=v;
    edge[op].f=f;
    edge[op].next=map[u].next;
    map[u].next=op++;
}
bool bfs() {
    int head,tail,temp,mid,tt;
    head=0;
    tail=1;
    que[0]=1;
    memset(mark,0,sizeof(mark));
    mark[1]=1;
    while(head<tail) {
        temp=que[head++];
        if(temp==n) return 1;
        mid=map[temp].next;
        while(mid!=-1) {
            tt=edge[mid].p;
            if(edge[mid].f>0&&!mark[tt]) {
                mark[tt]=1;
                que[tail++]=tt;
                pre[tt]=mid;
            }
            mid=edge[mid].next;
        }
    }
    return 0;
}
int min(int a,int b) {
    return a<b?a:b;
}
void max_flow() {
    int temp,fin;
    fin=INF;
    temp=n;
    while(temp!=1) {
        fin=min(fin,edge[pre[temp]].f);
        temp=edge[pre[temp]^1].p;
    }
    temp=n;
    ans+=fin;
    while(temp!=1) {
        edge[pre[temp]].f-=fin;
        edge[pre[temp]^1].f+=fin;
        temp=edge[pre[temp]^1].p;
    }
}
int main() {
    int m,a,b,c,i;
    while(~scanf("%d%d",&m,&n)) {
        op=0;
        for(i=1; i<=n; i++)
            map[i].next=-1;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            add(a,b,c);
            add(b,a,0);
        }
        ans=0;
        while(bfs())
            max_flow();
        printf("%d\n",ans);
    }
    return 0;
}

Dinic
struct node {
    int p,f,next;
    node() {}
    node(int a,int b,int c):p(a),f(b),next(c) {}
} edge[410];
int map[210],que[210],dis[210],start[210];
bool mark[210];
int s,t,ans,op;
bool bfs() {
    int head=0,tail=1,i,tt,temp;
    for(i=s; i<=t; i++) dis[i]=-1;
    que[0]=s,dis[s]=0;
    while(head<tail) {
        temp=que[head++];
        for(i=map[temp]; ~i; i=edge[i].next) {
            tt=edge[i].p;
            if(dis[tt]==-1&&edge[i].f>0) dis[tt]=dis[temp]+1,que[tail++]=tt;
        }
    }
    return dis[t]!=-1;
}
int dfs(int pos,int f) {
    int i,temp,tt;
    mark[pos]=1;
    if(pos==t) {
        ans+=f;
        return f;
    }
    for(i=start[pos]; ~i; i=edge[i].next,start[pos]=i) {
        tt=edge[i].p;
        if(!mark[tt]&&dis[tt]==dis[pos]+1&&edge[i].f>0) {
            if(temp=dfs(tt,min(f,edge[i].f))) {
                edge[i].f-=temp;
                edge[i^1].f+=temp;
                return temp;
            }
        }
    }
    return 0;
}
void dinic() {
    ans=0;
    int i;
    while(bfs()) {
        for(i=s; i<=t; i++) start[i]=map[i];
        memset(mark,0,sizeof(mark));
        while(dfs(s,inf)) memset(mark,0,sizeof(mark));
    }
    printf("%d\n",ans);
}
void add(int a,int b,int c) {
    edge[op]=node(b,c,map[a]);
    map[a]=op++;
}
int main() {
    int n,m,i,a,b,c;
    while(~scanf("%d%d",&m,&n)) {
        for(i=1,op=0; i<=n; i++) map[i]=-1;
        s=1,t=n;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            add(a,b,c),add(b,a,0);
        }
        dinic();
    }
    return 0;
}
最小费用最大流
struct node {
    int p,next;
    int f,c;
    node() {}
    node(int a,int b,int cc,int d):p(a),f(b),c(cc),next(d) {}
} edge[100000],map[1010];
int op,s,t;
int dis[1010],pre[1010];
bool mark[1010];
int que[2000000];
int spfa() {
    int head,tail,temp,mid,i,tt;
    head=0;
    tail=1;
    que[0]=s;
    for(i=s; i<=t; i++) {
        dis[i]=inf;
        mark[i]=0;
    }
    dis[s]=0;
    mark[s]=1;
    while(head<tail) {
        temp=que[head++];
        mark[temp]=0;
        mid=map[temp].next;
        while(mid!=-1) {
            tt=edge[mid].p;
            if(dis[tt]>dis[temp]+edge[mid].c&&edge[mid].f>0) {
                dis[tt]=dis[temp]+edge[mid].c;
                pre[tt]=mid;
                if(!mark[tt]) {
                    mark[tt]=1;
                    que[tail++]=tt;
                }
            }
            mid=edge[mid].next;
        }
    }
    return dis[t];
}

void add(int a,int b,int f,int c) {
    edge[op]=node(b,f,c,map[a].next);
    map[a].next=op++;
}
void solve() {
    int ans,temp,cost,fin=inf;
    ans=0;
    while(1) {
        cost=spfa();
        if(cost==inf) break;
        temp=t;
        while(temp!=s) {
            fin=min(fin,edge[pre[temp]].f);
            temp=edge[pre[temp]^1].p;
        }
        ans+=cost;
        temp=t;
        while(temp!=s) {
            edge[pre[temp]].f-=fin;
            edge[pre[temp]^1].f+=fin;
            temp=edge[pre[temp]^1].p;
        }
    }
    printf("%d\n",ans);
}
int main() {
    int n,m,i,a,b,c;
    scanf("%d%d",&n,&m);
    op=0;
    s=0;
    t=n+1;
    for(i=s; i<=t; i++)
        map[i].next=-1;
    while(m--) {
        scanf("%d%d%d",&a,&b,&c);
        add(a,b,1,c);
        add(b,a,0,-c);
        add(b,a,1,c);
        add(a,b,0,-c);
    }
    add(s,1,2,0);
    add(1,s,0,0);
    add(n,t,2,0);
    add(t,n,0,0);
    solve();
}
有上下界的最大流EK
const int inf=0x7FFFFFFF,v_size=200,e_size=1000;
struct node {
    int p,next,f,c;
    node() {}
    node(int a,int b,int cc,int d):p(a),f(b),c(cc),next(d) {}
} edge[e_size],map[v_size];
int op,ans,ss,tt;
void add(int a,int b,int c,int d) {
    edge[op]=node(b,c,d,map[a].next);
    map[a].next=op++;
}
int que[1000000],pre[v_size],st[v_size],se[v_size];
bool mark[v_size];
bool bfs() {
    int head,tail,temp,i,mid;
    head=0;
    tail=1;
    que[0]=ss;
    memset(mark,0,sizeof(mark));
    mark[ss]=1;
    while(head<tail) {
        temp=que[head++];
        if(temp==tt) return 1;
        for(i=map[temp].next; i+1; i=edge[i].next) {
            mid=edge[i].p;
            if(!mark[mid]&&edge[i].f>0) {
                pre[mid]=i;
                mark[mid]=1;
                que[tail++]=mid;
            }
        }
    }
    return 0;
}
void solve() {
    int temp,fin;
    ans=0;
    while(bfs()) {
        fin=inf;
        temp=tt;
        while(temp!=ss) {
            fin=min(fin,edge[pre[temp]].f);
            temp=edge[pre[temp]^1].p;
        }
        ans+=fin;
        temp=tt;
        while(temp!=ss) {
            edge[pre[temp]].f-=fin;
            edge[pre[temp]^1].f+=fin;
            temp=edge[pre[temp]^1].p;
        }
    }
}
int main() {
    int n,m,s,t,i,flow,b,a,c,d,j;
    while(scanf("%d%d",&n,&m)!=EOF) {
        s=1;
        t=n;
        ss=0;
        tt=n+1;
        for(i=ss; i<=tt; i++)
            map[i].next=-1;
        op=0;
        flow=0;
        memset(st,0,sizeof(st));
        memset(se,0,sizeof(se));
        while(m--) {
            scanf("%d%d%d%d",&a,&b,&c,&d);
            st[a]+=c;
            se[b]+=c;
            flow+=c;
            add(a,b,d-c,d-c);
            add(b,a,0,0);
        }
        for(i=1; i<=n; i++) {
            add(ss,i,se[i],se[i]);
            add(i,ss,0,0);
            add(i,tt,st[i],st[i]);
            add(tt,i,0,0);
        }
        add(t,s,inf,inf);
        add(s,t,0,0);
        solve();
        if(ans!=flow)
            printf("no answer\n");
        else {
            map[s].next=edge[map[s].next].next;
            map[t].next=edge[map[t].next].next;
            for(i=s; i<=t; i++)
                for(j=1; j<=2; j++)
                    map[i].next=edge[map[i].next].next;
            ss=s;
            tt=t;
            solve();
            ans=0;
            for(i=map[s].next; i+1; i=edge[i].next)
                ans+=edge[i].c-edge[i].f;
            printf("%d\n",ans+st[s]);
        }
    }
}
有上下界的最小费最大流
#include<iostream>
using namespace std;
const int inf=0x7FFFFFFF;
const int v_size=110,e_size=1000;
struct node {
    int p,next,f,c;
    node() {}
    node(int a,int b,int cc,int d):p(a),f(b),c(cc),next(d) {}
} edge[e_size],map[v_size];
int dis[v_size],pre[v_size];
bool mark[v_size];
int op,ss,tt,ans,n;
int que[1000000];
void add(int a,int b,int c,int d) {
    edge[op]=node(b,c,d,map[a].next);
    map[a].next=op++;
}
int spfa() {
    int head,tail,i,mid,temp;
    head=0;
    tail=1;
    que[0]=ss;
    for(i=0; i<=n+3; i++) {
        dis[i]=inf;
        mark[i]=0;
    }
    mark[ss]=1;
    dis[ss]=0;
    while(head<tail) {
        temp=que[head++];
        for(i=map[temp].next; i+1; i=edge[i].next) {
            if(edge[i].f<=0) continue;
            mid=edge[i].p;
            if(dis[mid]>dis[temp]+edge[i].c) {
                dis[mid]=dis[temp]+edge[i].c;
                pre[mid]=i;
                if(!mark[mid]) {
                    mark[mid]=1;
                    que[tail++]=mid;
                }
            }
        }
        mark[temp]=0;
    }
    return dis[tt];
}
void solve() {
    int fin,cost,temp;
    while(1) {
        cost=spfa();
        if(cost==inf) break;
        temp=tt;
        fin=inf;
        while(temp!=ss) {
            fin=min(fin,edge[pre[temp]].f);
            temp=edge[pre[temp]^1].p;
        }
        temp=tt;
        ans+=fin*cost;
        while(temp!=ss) {
            edge[pre[temp]].f-=fin;
            edge[pre[temp]^1].f+=fin;
            temp=edge[pre[temp]^1].p;
        }
    }
}
int main() {
    int t,s,a,sum,m,i,b,c;
    while(scanf("%d",&n)!=EOF) {
        s=n;
        t=n+1;
        ss=n+2;
        tt=n+3;
        for(i=0; i<=n+3; i++)
            map[i].next=-1;
        op=0;
        sum=0;
        m=n-1;
        for(i=0; i<n; i++) {
            scanf("%d",&a);
            sum+=a;
            add(s,i,a,0);
            add(i,s,0,0);
        }
        sum/=n;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            add(a,b,inf,c);
            add(b,a,0,-c);
            add(b,a,inf,c);
            add(a,b,0,-c);
        }
        for(i=0; i<n; i++) {
            add(i,t,1,0);
            add(t,i,0,0);
            add(i,tt,sum,0);
            add(tt,i,0,0);
        }
        add(ss,t,sum*n,0);
        add(t,ss,0,0);
        add(t,s,inf,0);
        add(s,t,0,0);
        ans=0;
        solve();
        map[s].next=edge[map[s].next].next;
        map[t].next=edge[map[t].next].next;
        map[t].next=edge[map[t].next].next;
        for(i=0; i<n; i++)
            map[i].next=edge[map[i].next].next;
        tt=t;
        ss=s;
        solve();
        printf("%d\n",ans);
    }
}
//有上下界的最小流EK版
const int inf=0x7FFFFFFF;
struct node {
    int p,next,f;
    node() {}
    node(int a,int b,int c):p(a),f(b),next(c) {}
} edge[10000],map[60];
bool mark[60];
int que[60],pre[60],st[60],se[60];
int n,op,ss,tt,ans,sum;
void add(int a,int b,int c) {
    edge[op]=node(b,c,map[a].next);
    map[a].next=op++;
}
bool bfs() {
    int head,tail,temp,mid,i;
    memset(mark,0,sizeof(mark));
    head=0;
    tail=1;
    que[0]=ss;
    mark[ss]=1;
    while(head<tail) {
        temp=que[head++];
        if(temp==tt) return 1;
        for(i=map[temp].next; i+1; i=edge[i].next) {
            mid=edge[i].p;
            if(!mark[mid]&&edge[i].f>0) {
                pre[mid]=i;
                mark[mid]=1;
                que[tail++]=mid;
            }
        }
    }
    return 0;
}
void solve() {
    int fin,temp;
    ans=0;
    while(bfs()) {
        temp=tt;
        fin=inf;
        while(temp!=ss) {
            fin=min(fin,edge[pre[temp]].f);
            temp=edge[pre[temp]^1].p;
        }
        temp=tt;
        ans+=fin;
        while(temp!=ss) {
            edge[pre[temp]].f-=fin;
            edge[pre[temp]^1].f+=fin;
            temp=edge[pre[temp]^1].p;
        }
    }
    if(ans!=sum) printf("impossible\n");
    else printf("%d\n",inf-edge[op-2].f);
}
int main() {
    int m,n,s,t,i,a,b,c;
    char stra[10],strb[10];
    while(scanf("%d%d",&n,&m),n+m) {
        s=0;
        t=n+1;
        ss=n+2;
        tt=n+3;
        for(i=0; i<=n+3; i++)
            map[i].next=-1;
        op=0;
        sum=0;
        memset(st,0,sizeof(st));
        memset(se,0,sizeof(se));
        while(m--) {
            scanf("%s%s%d",stra,strb,&c);
            sum+=c;
            if(stra[0]=='+') a=0;
            if(strb[0]=='-') b=n+1;
            if(isdigit(stra[0])) a=atoi(stra);
            if(isdigit(strb[0])) b=atoi(strb);
            add(a,b,inf);
            add(b,a,0);
            st[a]+=c;
            se[b]+=c;
        }
        for(i=0; i<=n+1; i++) {
            add(ss,i,se[i]);
            add(i,ss,0);
            add(i,tt,st[i]);
            add(tt,i,0);
        }
        add(t,s,inf);
        add(s,t,0);
        solve();
    }
}
KM最大权二分图匹配
int match[110],xl[110],yl[110],slack[110],dis[110][110];
bool vx[110],vy[110];
struct node {
    int x,y;
} hp[110],mp[110];
int hn,mn;
bool path(int pos) {
    int i,temp;
    vx[pos]=1;
    for(i=1; i<=hn; i++) {
        if(!vy[i]) {
            temp=xl[pos]+yl[i]-dis[pos][i];
            if(temp==0) {
                vy[i]=1;
                if(match[i]==-1||path(match[i])) {
                    match[i]=pos;
                    return 1;
                }
            } else slack[i]=min(slack[i],temp);
        }
    }
    return 0;
}
void km() {
    int i,j,ans,fin;
    for(i=1; i<=mn; i++) {
        xl[i]=-inf;
        for(j=1; j<=hn; j++)
            xl[i]=max(xl[i],dis[i][j]);
    }
    for(i=1; i<=hn; i++)
        yl[i]=0;
    memset(match,-1,sizeof(match));
    for(i=1; i<=mn; i++) {
        for(j=1; j<=hn; j++)
            slack[j]=inf;
        while(1) {
            memset(vx,0,sizeof(vx));
            memset(vy,0,sizeof(vy));
            if(path(i)) break;
            fin=inf;
            for(j=1; j<=hn; j++)
                if(!vy[j]) fin=min(fin,slack[j]);
            for(j=1; j<=mn; j++)
                if(vx[j]) xl[j]-=fin;
            for(j=1; j<=hn; j++) {
                if(vy[j]) yl[j]+=fin;
                else slack[j]-=fin;
            }
        }
    }
    ans=0;
    for(i=1; i<=hn; i++)
        if(match[i]!=-1) ans+=(0-dis[match[i]][i]);
    printf("%d\n",ans);
}
int main() {
    int n,m,i,j;
    char str[110];
    while(scanf("%d%d",&n,&m),n+m) {
        mn=0;
        hn=0;
        for(i=0; i<n; i++) {
            scanf("%s",str);
            for(j=0; j<m; j++) {
                if(str[j]=='H') {
                    ++hn;
                    hp[hn].x=i;
                    hp[hn].y=j;
                }
                if(str[j]=='m') {
                    ++mn;
                    mp[mn].x=i;
                    mp[mn].y=j;
                }
            }
        }
        for(i=1; i<=mn; i++)
            for(j=1; j<=hn; j++)
                dis[i][j]=-(abs(mp[i].x-hp[j].x)+abs(mp[i].y-hp[j].y));
        km();
    }
}
差分约束
struct node {
    int s,len;
} rec[500010];
int p,head,tail,size;
int start[50010],next[500010],dis[50010];
int que[1000010];
bool mark[50010];
void spfa() {
    int mid,temp,i,tt;
    head=0;
    tail=1;
    que[0]=0;
    for(i=1; i<=size; i++)
        dis[i]=-inf;
    memset(mark,0,sizeof(mark));
    mark[0]=1;
    dis[0]=0;
    while(head<tail) {
        temp=que[head++];
        mid=start[temp];
        while(mid!=-1) {
            tt=rec[mid].s;
            if(dis[tt]<dis[temp]+rec[mid].len) {
                dis[tt]=dis[temp]+rec[mid].len;
                if(!mark[tt]) {
                    mark[tt]=1;
                    que[tail++]=tt;
                }
            }
            mid=next[mid];
        }
        mark[temp]=0;
    }
    printf("%d\n",dis[size]);
}
void ins(int a,int b,int c) {
    rec[p].s=b;
    rec[p].len=c;
    next[p]=start[a];
    start[a]=p;
    p++;
}
int main() {
    int n,a,b,c,i;
    while(~scanf("%d",&n)) {
        size=0;
        p=0;
        memset(start,-1,sizeof(start));
        while(n--) {
            scanf("%d%d%d",&a,&b,&c);
            size=max(size,b+1);
            ins(a,b+1,c);
        }
        for(i=0; i<size; i++) {
            ins(i,i+1,0);
            ins(i+1,i,-1);
        }
        spfa();
    }
}
强连通分量tarjan算法
const int inf=0x7FFFFFFF;
struct node {
    int p,next;
    node() {}
    node(int a,int b):p(a),next(b) {}
} edge[1010];
int map[210],low[210],part[210],stack[210];
int top,op,ans,n;
void dfs(int pos) {
    int cur=low[pos]=op++,i,tt;
    stack[++top]=pos;
    for(i=map[pos]; i+1; i=edge[i].next) {
        tt=edge[i].p;
        if(low[tt]==-1) dfs(tt);
        cur=min(cur,low[tt]);
    }
    if(cur<low[pos]) {
        low[pos]=cur;
        return ;
    }
    ans++;
    do {
        part[tt=stack[top--]]=ans;
        low[tt]=inf;
    } while(tt!=pos);
}
void add(int a,int b) {
    edge[op]=node(b,map[a]);
    map[a]=op++;
}
int main() {
    int m,a,b,i;
    while(scanf("%d%d",&n,&m),n+m) {
        for(i=1; i<=n; i++)
            map[i]=-1;
        op=0;
        while(m--) {
            scanf("%d%d",&a,&b);
            add(a,b);
        }
        for(i=1; i<=n; i++)
            low[i]=-1;
        ans=0;
        op=0;
        top=0;
        for(i=1; i<=n; i++)
            if(low[i]==-1) dfs(i);
        printf("%d\n",ans);
    }
}
RMQ求解LCA问题
struct node {
    int p,w,next;
    node() {}
    node(int a,int b,int c):p(a),w(b),next(c) {}
} edge[20010];
int p[20];
int dp[20010][20];
int map[10010],mark[10010],deep[20010],e[20010],r[10010],dis[10010];
int op,num;
void add(int a,int b,int l) {
    edge[op]=node(b,l,map[a]);
    map[a]=op++;
}
void dfs(int pos,int len,int cost) {
    int i,tt;
    mark[pos]=num;
    dis[pos]=cost;
    r[pos]=op;
    e[op]=pos;
    deep[op++]=len;
    for(i=map[pos]; i+1; i=edge[i].next) {
        tt=edge[i].p;
        if(mark[tt]) continue;
        dfs(tt,len+1,cost+edge[i].w);
        e[op]=pos;
        deep[op++]=len;
    }
}
int query(int x,int y) {
    int k=(int)(log((double)(y+1-x))/log(2.0));
    return deep[dp[x][k]]>deep[dp[y-p[k]+1][k]]?dp[y-p[k]+1][k]:dp[x][k];
}
void initrmq() {
    int j,i;
    for(i=0; i<op; i++) dp[i][0]=i;
    int k=(int)(log((double)op)/log(2.0));
    for(j=1; j<=k; j++)
        for(i=0; i<op; i++) {
            if(i+p[j]-1<op)
                dp[i][j]=deep[dp[i][j-1]]>deep[dp[i+p[j-1]][j-1]]?dp[i+p[j-1]][j-1]:dp[i][j-1];
            else break;
        }
}
int main() {
    int n,m,c,a,b,l,i;
    for(i=1,p[0]=1; i<=14; i++)
        p[i]=p[i-1]<<1;
    while(~scanf("%d%d%d",&n,&m,&c)) {
        op=0;
        for(i=1; i<=n; i++) {
            mark[i]=0;
            map[i]=-1;
        }
        while(m--) {
            scanf("%d%d%d",&a,&b,&l);
            add(a,b,l);
            add(b,a,l);
        }
        op=0;
        num=0;
        for(i=1; i<=n; i++) {
            num++;
            if(!mark[i]) dfs(i,0,0);
        }
        initrmq();
        while(c--) {
            scanf("%d%d",&a,&b);
            if(mark[a]!=mark[b]) printf("Not connected\n");
            else printf("%d\n",dis[a]+dis[b]-2*dis[e[query(min(r[a],r[b]),max(r[a],r[b]))]]);
        }
    }
}
最小树形图
const double inf=1e10;
struct node {
    double x,y;
    node() {}
    node(double a,double b):x(a),y(b) {}
} rec[110];
bool mark[110],vis[110];
int pre[110];
int n;
double map[110][110];
double dis(node a,node b) {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
void dfs(int pos) {
    int i;
    mark[pos]=1;
    for(i=1; i<=n; i++) {
        if(mark[i]||map[pos][i]==inf) continue;
        dfs(i);
    }
}
bool can() {
    int i,judge=0;
    memset(mark,0,sizeof(mark));
    for(i=1; i<=n; i++)
        if(!mark[i]) {
            judge++;
            dfs(i);
        }
    return judge==1;
}
void solve() {
    double ans;
    int i,j,k;
    memset(mark,0,sizeof(mark));
    while(1) {
        for(i=2; i<=n; i++) {
            if(mark[i]) continue;
            pre[i]=i;
            map[i][i]=inf;
            for(j=1; j<=n; j++) {
                if(mark[j]) continue;
                if(map[j][i]<map[pre[i]][i]) pre[i]=j;
            }
        }
        for(i=2; i<=n; i++) {
            if(mark[i]) continue;
            memset(vis,0,sizeof(vis));
            j=i;
            while(!vis[j]&&j!=1) {
                vis[j]=1;
                j=pre[j];
            }
            if(j==1) continue;
            i=j;
            ans+=map[pre[i]][i];
            for(j=pre[i]; j!=i; j=pre[j]) {
                ans+=map[pre[j]][j];
                mark[j]=1;
            }
            for(j=1; j<=n; j++) {
                if(mark[j]) continue;
                if(map[j][i]!=inf) map[j][i]-=map[pre[i]][i];
            }
            for(j=pre[i]; j!=i; j=pre[j]) {
                for(k=1; k<=n; k++) {
                    if(mark[k]) continue;
                    map[i][k]=min(map[i][k],map[j][k]);
                    if(map[k][j]!=inf)
                        map[k][i]=min(map[k][i],map[k][j]-map[pre[j]][j]);
                }
            }
            break;
        }
        if(i>n) {
            for(j=2; j<=n; j++)
                if(!mark[j]) ans+=map[pre[j]][j];
            break;
        }
    }
    printf("%.2lf\n",ans);
}
int main() {
    int i,m,a,b,j;
    while(~scanf("%d%d",&n,&m)) {
        for(i=1; i<=n; i++)
            for(j=1; j<=n; j++)
                map[i][j]=inf;
        for(i=1; i<=n; i++)
            scanf("%lf%lf",&rec[i].x,&rec[i].y);
        while(m--) {
            scanf("%d%d",&a,&b);
            map[a][b]=dis(rec[a],rec[b]);
        }
        if(!can()) {
            printf("poor snoopy\n");
            continue;
        }
        solve();
    }
}
Floyd求最小环
int main() {
    int i,j,k,ans,a,b,c,n,m;
    while(scanf("%d%d",&n,&m)==2) {
        for(i=1; i<=n; i++)
            for(j=1; j<=n; j++)
                if(i==j) lin[i][j]=0;
                else lin[i][j]=inf;
        while(m--) {
            scanf("%d%d%d",&a,&b,&c);
            if(lin[a][b]>c) lin[a][b]=lin[b][a]=c;
        }
        ans=inf;
        for(i=1; i<=n; i++)
            for(j=1; j<=n; j++)
                dis[i][j]=lin[i][j];
        for(k=1; k<=n; k++) {
            for(i=1; i<=k-1; i++)
                for(j=i+1; j<=k-1; j++)
                    ans=min(ans,lin[i][j]+dis[i][k]+dis[k][j]);
            for(i=1; i<=n; i++)
                for(j=1; j<=n; j++)
                    lin[i][j]=min(lin[i][j],lin[i][k]+lin[k][j]);
        }
        if(ans==inf) printf("No solution.\n");
        else printf("%d\n",ans);
    }
}
最大权闭合子图
typedef long long ll;
const int MAXN=5010,MANM=200010;
const ll inf=1LL<<62;
struct node {
    int p,next;
    ll f;
    node() {}
    node(int a,ll b,int c):p(a),f(b),next(c) {}
} edge[MANM];
int map[MAXN],que[MAXN],dis[MAXN],start[MAXN];
bool mark[MAXN];
int s,t,op;
ll ans;
bool bfs() {
    int head=0,tail=1,i,tt,temp;
    for(i=s; i<=t; i++) dis[i]=-1;
    que[0]=s,dis[s]=0;
    while(head<tail) {
        temp=que[head++];
        for(i=map[temp]; ~i; i=edge[i].next) {
            tt=edge[i].p;
            if(dis[tt]==-1&&edge[i].f>0) dis[tt]=dis[temp]+1,que[tail++]=tt;
        }
    }
    return dis[t]!=-1;
}
ll dfs(int pos,ll f) {
    int i,tt;
    ll temp;
    mark[pos]=1;
    if(pos==t) {
        ans+=f;
        return f;
    }
    for(i=start[pos]; ~i; i=edge[i].next,start[pos]=i) {
        tt=edge[i].p;
        if(!mark[tt]&&dis[tt]==dis[pos]+1&&edge[i].f>0) {
            if(temp=dfs(tt,min(f,edge[i].f))) {
                edge[i].f-=temp;
                edge[i^1].f+=temp;
                return temp;
            }
        }
    }
    return 0;
}
ll sum;
int num;
void find(int pos) {
    num++;
    mark[pos]=1;
    int i;
    for(i=map[pos]; ~i; i=edge[i].next) {
        if(edge[i].f>0&&!mark[edge[i].p]) find(edge[i].p);
    }
}
void dinic() {
    ans=0;
    int i;
    while(bfs()) {
        for(i=s; i<=t; i++) start[i]=map[i],mark[i]=0;
        while(dfs(s,inf)) for(i=s; i<=t; i++) mark[i]=0;
    }
    num=0;
    for(i=s; i<=t; i++) mark[i]=0;
    find(s);
    printf("%d %lld\n",num-1,sum-ans);
}
void add(int a,int b,ll c) {
    edge[op]=node(b,c,map[a]);
    map[a]=op++;
}
int main() {
    int n,m,i,a,b;
    ll w;
    while(~scanf("%d%d",&n,&m)) {
        s=0,t=n+1;
        for(i=s,op=0; i<=t; i++) map[i]=-1;
        sum=0;
        for(i=1; i<=n; i++) {
            scanf("%lld",&w);
            if(w>0) sum+=w,add(s,i,w),add(i,s,0);
            else add(i,t,(0-w)),add(t,i,0);
        }
        while(m--) {
            scanf("%d%d",&a,&b);
            add(a,b,inf);
            add(b,a,0);
        }
        dinic();
    }
    return 0;
}
最大密度子图 二分法+dinic
const double eps=1e-6,inf=1e10;
struct node {
    int p,next;
    double w;
    node() {}
    node(int a,double b,int c):p(a),w(b),next(c) {}
} edge[10000];
struct gg {
    int ss,tt;
} rec[1010];
int dis[110],map[110],que[110],start[110],d[110];
bool mark[110];
int n,m,s,t,op;
bool bfs() {
    int head=0,tail=1,i,temp,tt;
    for(i=s; i<=t; i++) dis[i]=-1;
    dis[s]=0,que[0]=s;
    while(head<tail) {
        temp=que[head++];
        for(i=map[temp]; ~i; i=edge[i].next) {
            tt=edge[i].p;
            if(dis[tt]==-1&&edge[i].w>0) dis[tt]=dis[temp]+1,que[tail++]=tt;
        }
    }
    return dis[t]!=-1;
}
double dfs(int pos,double f,double &ans) {
    int i,tt;
    double temp;
    mark[pos]=1;
    if(pos==t) {
        ans+=f;
        return f;
    }
    for(i=start[pos]; ~i; i=edge[i].next,start[pos]=i) {
        tt=edge[i].p;
        if(!mark[tt]&&dis[tt]==dis[pos]+1&&edge[i].w>0) {
            temp=dfs(tt,min(f,edge[i].w),ans);
            if(temp>0) {
                edge[i].w-=temp,edge[i^1].w+=temp;
                return temp;
            }
        }
    }
    return 0;
}
double dinic() {
    double ans=0;
    int i;
    while(bfs()) {
        for(i=s; i<=t; i++) start[i]=map[i];
        memset(mark,0,sizeof(mark));
        while(dfs(s,inf,ans)>eps) memset(mark,0,sizeof(mark));
    }
    return ans;
}
void add(int a,int b,double c) {
    edge[op]=node(b,c,map[a]);
    map[a]=op++;
}
bool check(double g) {
    int i;
    for(i=s,op=0; i<=t; i++) map[i]=-1;
    for(i=1; i<=n; i++) add(s,i,m),add(i,s,0),add(i,t,(double)m+2*g-(double)d[i]),add(t,i,0);
    for(i=1; i<=m; i++) add(rec[i].ss,rec[i].tt,1),add(rec[i].tt,rec[i].ss,1);
    return (double)m*(double)n>dinic()+eps;
}
int main() {
    int i;
    double mid,low,high;
    scanf("%d%d",&n,&m);
    if(m==0) {
        printf("0.0000\n");
        return 0;
    }
    for(i=1; i<=n; i++) d[i]=0;
    for(i=1; i<=m; i++) scanf("%d%d",&rec[i].ss,&rec[i].tt),d[rec[i].ss]++,d[rec[i].tt]++;
    low=1.0/(double)n,high=(double)m;
    s=0,t=n+1;
    while(high-low>eps) {
        mid=(high+low)/2;
        if(check(mid)) low=mid+eps;
        else high=mid-eps;
    }
    printf("%.4lf\n",mid);
    return 0;
}
混合图欧拉回路 hdu 3472
const int inf=0x3fffffff,Esize=10000,Vsize=30,s=26,t=27;
struct node {
    int p,f,next;
    node() {}
    node(int a,int b,int c):p(a),f(b),next(c) {}
} edge[Esize];
int op;
int start[Vsize],map[Vsize],dis[Vsize],indegree[Vsize],outdegree[Vsize],que[Vsize];
bool lin[Vsize][Vsize],in[Vsize],mark[Vsize];
void search(int pos) {
    mark[pos]=1;
    int i;
    for(i=0; i<=25; i++) if(lin[pos][i]&&!mark[i]) search(i);
}
bool check() {
    int i,judge=0;
    for(i=0; i<=25; i++) {
        if(in[i]&&!mark[i]) {
            judge++;
            if(judge>1) return 0;
            search(i);
        }
    }
    return 1;
}
void add(int a,int b,int c) {
    edge[op]=node(b,c,map[a]);
    map[a]=op++;
}
bool bfs() {
    int i,temp,tt,head=0,tail=1;
    for(i=0; i<=t; i++) dis[i]=-1;
    que[0]=s;
    dis[s]=0;
    while(head<tail) {
        temp=que[head++];
        for(i=map[temp]; ~i; i=edge[i].next) {
            tt=edge[i].p;
            if(dis[tt]==-1&&edge[i].f>0) dis[tt]=dis[temp]+1,que[tail++]=tt;
        }
    }
    return dis[t]!=-1;
}
int dfs(int pos,int flow,int &ans) {
    int tt,temp,i;
    mark[pos]=1;
    if(pos==t) {
        ans+=flow;
        return flow;
    }
    for(i=start[pos]; ~i; i=edge[i].next,start[pos]=i) {
        tt=edge[i].p;
        if(!mark[tt]&&dis[tt]==dis[pos]+1&&edge[i].f>0) {
            temp=dfs(tt,min(flow,edge[i].f),ans);
            if(temp>0) {
                edge[i].f-=temp;
                edge[i^1].f+=temp;
                return temp;
            }
        }
    }
    return 0;
}
int max_flow() {
    int ans=0,i;
    while(bfs()) {
        memset(mark,0,sizeof(mark));
        for(i=0; i<=t; i++) start[i]=map[i];
        while(dfs(s,inf,ans)) memset(mark,0,sizeof(mark));
    }
    return ans;
}
int main() {
    bool judge;
    int T,n,type,a,b,ss,cases=0,tt,i,flow,val,sum;
    char str[100];
    scanf("%d",&T);
    while(T--) {
        scanf("%d",&n);
        op=0;
        memset(map,-1,sizeof(map));
        memset(lin,0,sizeof(lin));
        memset(in,0,sizeof(in));
        memset(indegree,0,sizeof(indegree));
        memset(outdegree,0,sizeof(outdegree));
        while(n--) {
            scanf("%s%d",str,&type);
            a=str[0]-'a';
            b=str[strlen(str)-1]-'a';
            if(type) add(a,b,1),add(b,a,0);
            lin[a][b]=lin[b][a]=1;
            in[a]=in[b]=1;
            outdegree[a]++,indegree[b]++;
        }
        printf("Case %d: ",++cases);
        if(!check()) {
            printf("Poor boy!\n");
            continue;
        }
        ss=-1,tt=-1;
        judge=1;
        sum=0;
        for(i=0; i<=25; i++) {
            if(!in[i]) continue;
            val=indegree[i]-outdegree[i];
            if(val%2) {
                if(val<0) {
                    if(ss==-1) ss=i;
                    else judge=0;
                } else {
                    if(tt==-1) tt=i;
                    else judge=0;
                }
            } else {
                if(val<0) {
                    add(s,i,-val/2);
                    add(i,s,0);
                    sum+=(-val/2);
                } else if(val>0) {
                    add(i,t,val/2);
                    add(t,i,0);
                }
            }
        }
        if((ss==-1&&tt!=-1)||(ss!=-1&&tt==-1)) judge=0;
        if(!judge) {
            printf("Poor boy!\n");
            continue;
        }
        if(ss!=-1&&tt!=-1) {
            add(ss,tt,1);
            add(tt,ss,0);
            outdegree[ss]++,indegree[tt]++;
            val=indegree[ss]-outdegree[ss];
            if(val<0) {
                add(s,ss,-val/2);
                add(ss,s,0);
                sum+=(-val/2);
            } else if(val>0) {
                add(ss,t,val/2);
                add(t,ss,0);
            }
            val=indegree[tt]-outdegree[tt];
            if(val<0) {
                add(s,tt,-val/2);
                add(tt,s,0);
                sum+=(-val/2);
            } else if(val>0) {
                add(tt,t,val/2);
                add(t,tt,0);
            }
        }
        flow=max_flow();
        if(flow==sum) printf("Well done!\n");
        else printf("Poor boy!\n");
    }
    return 0;
}
