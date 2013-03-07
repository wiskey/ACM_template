//hash模板
unsigned int SDBMHash(char *str) {
    unsigned int hash = 0;
    while (*str) {
        hash = (*str++) + (hash << 6) + (hash << 16) - hash;
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int RSHash(char *str) {
    unsigned int b = 378551;
    unsigned int a = 63689;
    unsigned int hash = 0;
    while (*str) {
        hash = hash * a + (*str++);
        a *= b;
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int JSHash(char *str) {
    unsigned int hash = 1315423911;
    while (*str) {
        hash ^= ((hash << 5) + (*str++) + (hash >> 2));
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int PJWHash(char *str) {
    unsigned int BitsInUnignedInt = (unsigned int)(sizeof(unsigned int) * 8);
    unsigned int ThreeQuarters	= (unsigned int)((BitsInUnignedInt  * 3) / 4);
    unsigned int OneEighth		= (unsigned int)(BitsInUnignedInt / 8);
    unsigned int HighBits		 = (unsigned int)(0xFFFFFFFF) << (BitsInUnignedInt - OneEighth);
    unsigned int hash			 = 0;
    unsigned int test			 = 0;
    while (*str) {
        hash = (hash << OneEighth) + (*str++);
        if ((test = hash & HighBits) != 0) {
            hash = ((hash ^ (test >> ThreeQuarters)) & (~HighBits));
        }
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int ELFHash(char *str) {
    unsigned int hash = 0;
    unsigned int x	= 0;
    while (*str) {
        hash = (hash << 4) + (*str++);
        if ((x = hash & 0xF0000000L) != 0) {
            hash ^= (x >> 24);
            hash &= ~x;
        }
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int BKDRHash(char *str) {
    unsigned int seed = 131;
    unsigned int hash = 0;
    while (*str) {
        hash = hash * seed + (*str++);
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int DJBHash(char *str) {
    unsigned int hash = 5381;
    while (*str) {
        hash += (hash << 5) + (*str++);
    }
    return (hash & 0x7FFFFFFF);
}
unsigned int APHash(char *str) {
    unsigned int hash = 0;
    int i;
    for (i=0; *str; i++) {
        if ((i & 1) == 0) {
            hash ^= ((hash << 7) ^ (*str++) ^ (hash >> 3));
        } else {
            hash ^= (~((hash << 11) ^ (*str++) ^ (hash >> 5)));
        }
    }
    return (hash & 0x7FFFFFFF);
}
//KMP模板
char W[10010],T[1000010];
int next[10010];
int wlen,tlen;
inline void kmp() {
    int i,k,j,ans;
    next[0]=-1;
    i=0;
    k=-1;
    while(i<wlen) {
        if(k==-1||W[k]==W[i]) next[++i]=++k;
        else k=next[k];
    }
    i=-1;
    j=-1;
    ans=0;
    while(i<tlen) {
        if(j==-1||W[j]==T[i]) {
            i++;
            j++;
        } else j=next[j];
        if(j==wlen) ans++;
    }
    printf("%d\n",ans);
}
int main() {
    int t;
    scanf("%d",&t);
    while(t--) {
        scanf("%s%s",W,T);
        wlen=strlen(W);
        tlen=strlen(T);
        kmp();
    }
}
//AC自动机
struct node {
    node *fail;
    node *next[26];
    int count;
    node() {
        fail=NULL;
        memset(next,NULL,sizeof(next));
        count=0;
    }
}*que[1000000],*root;
char str[100010];
void insert(char *s) {
    node *p=root;
    int num,i=0;
    while(s[i]) {
        num=s[i]-'a';
        if(p->next[num]==NULL) p->next[num]=new node();
        p=p->next[num];
        i++;
    }
    p->count++;
}
void build() {
    int i,head=0,tail=1;
    node *p,*temp;
    root->fail=NULL;
    que[0]=root;
    while(head<tail) {
        temp=que[head++];
        for(i=0; i<26; i++) {
            if(temp->next[i]!=NULL) {
                if(temp==root) temp->next[i]->fail=root;
                else {
                    p=temp->fail;
                    while(p!=NULL) {
                        if(p->next[i]!=NULL) {
                            temp->next[i]->fail=p->next[i];
                            break;
                        }
                        p=p->fail;
                    }
                    if(p==NULL) temp->next[i]->fail=root;
                }
                que[tail++]=temp->next[i];
            }
        }
    }
}
void search() {
    int i=0,cnt=0,num;
    node *p=root,*temp;
    while(str[i]) {
        num=str[i]-'a';
        while(p->next[num]==NULL&&p!=root) p=p->fail;
        p=p->next[num];
        if(p==NULL) p=root;
        temp=p;
        while(temp!=root&&temp->count!=-1) {
            cnt+=temp->count;
            temp->count=-1;
            temp=temp->fail;
        }
        i++;
    }
    printf("%d\n",cnt);
}
int main() {
    int t,n;
    char keyword[60];
    scanf("%d",&t);
    while(t--) {
        scanf("%d",&n);
        root=new node();
        while(n--) {
            scanf("%s",keyword);
            insert(keyword);
        }
        build();
        scanf("%s",str);
        search();
    }
}
//后缀数组
int rank[MAX],trank[MAX],f[MAX],sa[MAX],tsa[MAX],cnt[MAX],h[MAX],dp[MAX][13],p[13];
int n;
char str[MAX];
bool equal(int tmp[],int a,int b,int k) {
    if(a+k>=n&&b+k>=n) return tmp[a]==tmp[b];
    if(a+k<n&&b+k<n) return tmp[a]==tmp[b]&&tmp[a+k]==tmp[b+k];
    return 0;
}
void creat_suffix_array() {
    int i,k;
    memset(cnt,0,sizeof(cnt));
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
void init_rmq() {
    int i,j;
    for(i=0; i<n; i++) dp[i][0]=h[i];
    int t=(int)(log((double)(n-1))/log(2.0));
    for(i=1; i<=t; i++)
        for(j=0; j<n; j++) {
            if(j+(1<<i)>=n) break;
            dp[j][i]=min(dp[j][i-1],dp[j+(1<<(i-1))][i-1]);
        }
}
int query(int a,int b) {
    int ta=rank[a],tb=rank[b];
    if(ta>tb) ta^=tb,tb^=ta,ta^=tb;
    int t=(int)(log((double)(tb-ta))/log(2.0));
    return min(dp[ta][t],dp[tb-(1<<t)][t]);
}
int main() {
    int len,i,ans,l,r,k,r1,r2,pos;
    for(p[0]=1,i=1; i<=11; i++) p[i]=p[i-1]<<1;
    while(~scanf("%s",str)) {
        len=strlen(str);
        str[len]='#';
        n=len+len+1;
        for(i=len+1; i<n; i++) str[i]=str[n-1-i];
        str[n]=0;
        creat_suffix_array();
        creat_height();
        init_rmq();
        ans=0;
        for(i=0; i<len; i++) {
            if(i) {
                k=query(i,n-i);
                if(ans<2*k) ans=2*k,pos=i-k;
            }
            k=query(i,n-i-1);
            if(ans<2*k-1) ans=2*k-1,pos=i-k+1;
        }
        for(i=pos; i<pos+ans; i++) printf("%c",str[i]);
        printf("\n");
    }
    return 0;
}
//堆模板
struct node {
    int id,w;
    char str[12];
};
class heap {
private:
    node que[size];
    int max;
    int pos[size];
public:
    void reset() {
        max=0;
    }
    node get(int num) {
        return que[pos[num]];
    }
    bool check(node temp,node mid) {
        return (temp.w>mid.w)||(temp.w==mid.w&&temp.id<mid.id);
    }
    void push(node temp) {
        int i;
        if(max==0) {
            que[0]=temp;
            pos[temp.id]=0;
        } else {
            for(i=max; i>0; i=(i-1)/2) {
                if(check(que[(i-1)/2],temp)) break;
                que[i]=que[(i-1)/2];
                pos[que[i].id]=i;
            }
            que[i]=temp;
            pos[que[i].id]=i;
        }
        max++;
    }
    node top() {
        return que[0];
    }
    void dele(int num) {
        int temp=pos[num];
        node mid;
        max--;
        que[temp]=que[max];
        pos[que[temp].id]=temp;
        mid=que[temp];
        int i,child;
        i=temp;
        child=2*i+1;
        while(child<max) {
            if(child<max-1&&check(que[child+1],que[child])) child++;
            if(check(mid,que[child])) break;
            else {
                que[i]=que[child];
                pos[que[i].id]=i;
                i=child;
                child=2*i+1;
            }
        }
        que[i]=mid;
        pos[que[i].id]=i;
    }
    void change(int num,int weight) {
        int temp=pos[num];
        if(weight>que[temp].w) {
            que[temp].w=weight;
            up(temp);
        } else if(weight<que[temp].w) {
            que[temp].w=weight;
            low(temp);
        }
    }
    void up(int num) {
        int i;
        node temp;
        temp=que[num];
        for(i=num; i>0; i=(i-1)/2) {
            if(check(que[(i-1)/2],temp)) break;
            que[i]=que[(i-1)/2];
            pos[que[i].id]=i;
        }
        que[i]=temp;
        pos[que[i].id]=i;
    }
    void low(int num) {
        node temp;
        temp=que[num];
        int child,i;
        i=num;
        child=2*i+1;
        while(child<max) {
            if(child<max-1&&check(que[child+1],que[child])) child++;
            if(check(temp,que[child])) break;
            else {
                que[i]=que[child];
                pos[que[i].id]=i;
                i=child;
                child=2*i+1;
            }
        }
        que[i]=temp;
        pos[que[i].id]=i;
    }
    bool empty() {
        return max==0;
    }
} H;
//线段树求逆序数
int T[4500000];
int ans;
void creat(int pos,int l,int r) {
    T[pos]=0;
    if(l==r) return;
    int mid=(l+r)>>1;
    creat(pos*2,l,mid);
    creat(pos*2+1,mid+1,r);
}
void insert(int pos,int l,int r,int a) {
    T[pos]++;
    if(l==r) return;
    int mid=(l+r)>>1;
    if(a<=mid) {
        ans+=T[pos*2+1];
        insert(pos*2,l,mid,a);
    } else insert(pos*2+1,mid+1,r,a);
}
int main() {
    int n,a,i;
    scanf("%d",&n);
    creat(1,1,n);
    ans=0;
    for(i=1; i<=n; i++) {
        scanf("%d",&a);
        insert(1,1,n,a);
    }
    printf("%d\n",ans);
}
//树状数组求逆序数
typedef long long ll;
const ll size=500010;
ll n;
struct node {
    ll num;
    ll rank;
} rec[size];
ll res[size],vis[size];
bool cmp(const node a,const node b) {
    return a.num<b.num;
}
ll lowbit(ll num) {
    return num&(-num);
}
void ins(ll num) {
    ll temp=num;
    while(temp>0) {
        vis[temp]++;
        temp-=lowbit(temp);
    }
}
ll sum(ll num) {
    ll temp=num;
    ll ans=0;
    while(temp<=n) {
        ans+=vis[temp];
        temp+=lowbit(temp);
    }
    return ans;
}
int main() {
    ll ans,i,j;
    while(scanf("%lld",&n),n) {
        for(i=0; i<=n-1; i++) {
            scanf("%lld",&rec[i].num);
            rec[i].rank=i;
        }
        sort(rec,rec+n,cmp);
        for(i=0; i<=n-1; i++) {
            rec[i].num=i+1;
        }
        for(i=0; i<=n-1; i++) {
            res[rec[i].rank]=rec[i].num;
        }
        ans=0;
        memset(vis,0,sizeof(vis));
        for(i=0; i<=n-1; i++) {
            ans+=sum(res[i]);
            ins(res[i]);
        }
        printf("%lld\n",ans);
    }
}
//线段树+扫描线+离散化求面积
struct node {
    int x,ly,hy,judge;
    node() {}
    node(int a,int b,int c,int d):x(a),ly(b),hy(c),judge(d) {}
    bool operator<(const node a)const {
        return x<a.x;
    }
} rec[2010];
struct gg {
    int y,rank;
    gg() {}
    gg(int a,int b):y(a),rank(b) {}
    bool operator<(const gg a)const {
        return y<a.y;
    }
} m[2010];
int map[2010];
struct tree {
    int cover,len;
    tree() {}
    tree(int a,int b):cover(a),len(b) {}
} T[10000];
void creat(int pos,int l,int r) {
    T[pos]=tree(0,0);
    if(l==r) return ;
    int mid=(l+r)>>1;
    creat(pos*2,l,mid);
    creat(pos*2+1,mid+1,r);
}
void deal(int pos,int l,int r) {
    if(T[pos].cover) {
        T[pos].len=map[r+1]-map[l];
        return ;
    }
    if(l==r) T[pos].len=0;
    else T[pos].len=T[pos*2].len+T[pos*2+1].len;
}
void insert(int pos,int l,int r,int start,int end,int oper) {
    if(l==start&&r==end) T[pos].cover+=oper;
    else {
        int mid=(l+r)>>1;
        if(end<=mid) insert(pos*2,l,mid,start,end,oper);
        else if(start>=mid+1) insert(pos*2+1,mid+1,r,start,end,oper);
        else {
            insert(pos*2,l,mid,start,mid,oper);
            insert(pos*2+1,mid+1,r,mid+1,end,oper);
        }
    }
    deal(pos,l,r);
}
int main() {
    int i,a,b,c,d,n,temp,j,mid,ans,lastlen;
    scanf("%d",&n);
    for(i=0; i<n; i++) {
        scanf("%d%d%d%d",&a,&b,&c,&d);
        rec[i<<1]=node(a,0,0,1);
        rec[(i<<1)^1]=node(c,0,0,-1);
        m[i<<1]=gg(b,i<<1);
        m[(i<<1)^1]=gg(d,(i<<1)^1);
    }
    sort(m,m+n*2);
    temp=-1;
    j=-1;
    for(i=0; i<n*2; i++) {
        if(temp==m[i].y) {
            mid=m[i].rank;
            if(mid&1) rec[mid].hy=rec[mid^1].hy=j;
            else rec[mid].ly=rec[mid^1].ly=j;
        } else {
            temp=m[i].y;
            mid=m[i].rank;
            map[++j]=temp;
            if(mid&1) rec[mid].hy=rec[mid^1].hy=j;
            else rec[mid].ly=rec[mid^1].ly=j;
        }
    }
    creat(1,0,j-1);
    ans=0;
    sort(rec,rec+2*n);
    for(i=0; i<n*2; i++) {
        insert(1,0,j-1,rec[i].ly,rec[i].hy-1,rec[i].judge);
        if(i) ans+=lastlen*(rec[i].x-rec[i-1].x);
        lastlen=T[1].len;
    }
    printf("%d\n",ans);
}
//线段树+扫描线+离散化求周长
struct node {
    int x,ly,hy,judge;
    node() {}
    node(int a,int b,int c,int d):x(a),ly(b),hy(c),judge(d) {}
    bool operator<(const node a)const {
        return x<a.x;
    }
} rec[10010];
struct gg {
    int y,rank;
    gg() {}
    gg(int a,int b):y(a),rank(b) {}
    bool operator<(const gg a)const {
        return y<a.y;
    }
} m[10010];
int map[10010];
struct tree {
    int cover,len,lis,ris,count;
    tree() {}
    tree(int a,int b,int c,int d,int e):cover(a),len(b),count(c),lis(d),ris(e) {}
} T[40010];
void creat(int pos,int l,int r) {
    T[pos]=tree(0,0,0,0,0);
    if(l==r) return ;
    int mid=(l+r)>>1;
    creat(pos*2,l,mid);
    creat(pos*2+1,mid+1,r);
}
void deal(int pos,int l,int r) {
    if(T[pos].cover) {
        T[pos].len=map[r+1]-map[l];
        T[pos].count=1;
        return ;
    }
    if(l==r) {
        T[pos].len=0;
        T[pos].count=0;
    } else {
        T[pos].len=T[pos*2].len+T[pos*2+1].len;
        T[pos].count=T[pos*2].count+T[pos*2+1].count;
        if(T[pos*2].ris&&T[pos*2+1].lis) T[pos].count--;
    }
}
void insert(int pos,int l,int r,int start,int end,int oper) {
    if(l==start) T[pos].lis+=oper;
    if(r==end) T[pos].ris+=oper;
    if(l==start&&r==end) T[pos].cover+=oper;
    else {
        int mid=(l+r)>>1;
        if(end<=mid) insert(pos*2,l,mid,start,end,oper);
        else if(start>=mid+1) insert(pos*2+1,mid+1,r,start,end,oper);
        else {
            insert(pos*2,l,mid,start,mid,oper);
            insert(pos*2+1,mid+1,r,mid+1,end,oper);
        }
    }
    deal(pos,l,r);
}
int main() {
    int i,a,b,c,d,n,temp,j,mid,ans,lastlen,lastnum;
    scanf("%d",&n);
    for(i=0; i<n; i++) {
        scanf("%d%d%d%d",&a,&b,&c,&d);
        rec[i<<1]=node(a,0,0,1);
        rec[(i<<1)^1]=node(c,0,0,-1);
        m[i<<1]=gg(b,i<<1);
        m[(i<<1)^1]=gg(d,(i<<1)^1);
    }
    sort(m,m+n*2);
    temp=-1;
    j=-1;
    for(i=0; i<n*2; i++) {
        if(temp==m[i].y) {
            mid=m[i].rank;
            if(mid&1) rec[mid].hy=rec[mid^1].hy=j;
            else rec[mid].ly=rec[mid^1].ly=j;
        } else {
            temp=m[i].y;
            mid=m[i].rank;
            map[++j]=temp;
            if(mid&1) rec[mid].hy=rec[mid^1].hy=j;
            else rec[mid].ly=rec[mid^1].ly=j;
        }
    }
    creat(1,0,j-1);
    ans=0;
    sort(rec,rec+2*n);
    lastlen=0;
    for(i=0; i<n*2; i++) {
        insert(1,0,j-1,rec[i].ly,rec[i].hy-1,rec[i].judge);
        ans+=abs(lastlen-T[1].len);
        lastlen=T[1].len;
        if(i) ans+=lastnum*(rec[i].x-rec[i-1].x)*2;
        lastnum=T[1].count;
    }
    printf("%d\n",ans);
}
//二维线段树求矩阵和
int n;
int Y1,Y2;
int T[2048][2048];
void creaty(int xx,int pos,int l,int r) {
    T[xx][pos]=0;
    if(l==r) return ;
    int mid=(l+r)>>1;
    creaty(xx,pos*2,l,mid);
    creaty(xx,pos*2+1,mid+1,r);
}
void creatx(int pos,int l,int r) {
    if(l==r) {
        creaty(pos,1,1,n);
        return;
    }
    int mid=(l+r)>>1;
    creaty(pos,1,1,n);
    creatx(pos*2,l,mid);
    creatx(pos*2+1,mid+1,r);
}
void inserty(int xx,int pos,int l,int r,int yy,int a) {
    T[xx][pos]+=a;
    if(l==r) return ;
    int mid=(l+r)>>1;
    if(yy<=mid) inserty(xx,pos*2,l,mid,yy,a);
    else inserty(xx,pos*2+1,mid+1,r,yy,a);
}
void insertx(int pos,int l,int r,int xx,int yy,int a) {
    if(l==r) {
        inserty(pos,1,1,n,yy,a);
        return ;
    }
    int mid=(l+r)>>1;
    inserty(pos,1,1,n,yy,a);
    if(xx<=mid) insertx(pos*2,l,mid,xx,yy,a);
    else insertx(pos*2+1,mid+1,r,xx,yy,a);
}
int sumy(int xx,int pos,int l,int r,int y1,int y2) {
    if(l==y1&&r==y2) return T[xx][pos];
    int mid=(l+r)>>1;
    if(y2<=mid) return sumy(xx,pos*2,l,mid,y1,y2);
    else if(y1>=mid+1) return sumy(xx,pos*2+1,mid+1,r,y1,y2);
    else return sumy(xx,pos*2,l,mid,y1,mid)+sumy(xx,pos*2+1,mid+1,r,mid+1,y2);
}
int sumx(int pos,int l,int r,int x1,int x2) {
    if(l==x1&&r==x2)
        return sumy(pos,1,1,n,Y1,Y2);
    int mid=(l+r)>>1;
    if(x2<=mid) return sumx(pos*2,l,mid,x1,x2);
    else if(x1>=mid+1) return sumx(pos*2+1,mid+1,r,x1,x2);
    else return sumx(pos*2,l,mid,x1,mid)+sumx(pos*2+1,mid+1,r,mid+1,x2);
}
int main() {
    int judge;
    int x1,x2,x,y,a;
    while(scanf("%d",&judge),judge-3) {
        if(judge==0) {
            scanf("%d",&n);
            creatx(1,1,n);
        }
        if(judge==1) {
            scanf("%d%d%d",&x,&y,&a);
            x++;
            y++;
            insertx(1,1,n,x,y,a);
        }
        if(judge==2) {
            scanf("%d%d%d%d",&x1,&Y1,&x2,&Y2);
            x1++;
            Y1++;
            x2++;
            Y2++;
            printf("%d\n",sumx(1,1,n,x1,x2));
        }
    }
}
//二维树状数组求矩阵和
const int inf=1999;
int rec[2000][2000];
int lowbit(int pos) {
    return pos&(-pos);
}
void insert(int x,int y,int a) {
    int tempx=x,tempy;
    while(tempx<inf) {
        tempy=y;
        while(tempy<inf) {
            rec[tempx][tempy]+=a;
            tempy+=lowbit(tempy);
        }
        tempx+=lowbit(tempx);
    }
}
int sum(int x,int y) {
    int tempx=x,tempy,ans=0;
    while(tempx>0) {
        tempy=y;
        while(tempy>0) {
            ans+=rec[tempx][tempy];
            tempy-=lowbit(tempy);
        }
        tempx-=lowbit(tempx);
    }
    return ans;
}
int main() {
    int t,n,j,i,x1,y1,x2,y2,a,ans;
    while((scanf("%d",&t),t)!=3) {
        if(t==0) {
            scanf("%d",&n);
            for(i=1; i<=n; i++)
                for(j=1; j<=n; j++)
                    rec[i][j]=0;
        }
        if(t==1) {
            scanf("%d%d%d",&x1,&y1,&a);
            x1++;
            y1++;
            insert(x1,y1,a);
        }
        if(t==2) {
            scanf("%d%d%d%d",&x1,&y1,&x2,&y2);
            x1++;
            y1++;
            x2++;
            y2++;
            ans=sum(x2,y2)-sum(x2,y1-1)-sum(x1-1,y2)+sum(x1-1,y1-1);
            printf("%d\n",ans);
        }
    }
}
//二维RMQ求矩阵最值
const int p[]= {1,1<<1,1<<2,1<<3,1<<4,1<<5,1<<6,1<<7,1<<8,1<<9};
int dp[9][9][310][310];
int n,m,r1,r2,c1,c2,ans;
int mylog(int temp) {
    int i;
    for(i=0;; i++)
        if(temp<p[i]) return i-1;
}
void redeal() {
    int logn=mylog(n);
    int logm=mylog(m);
    int i,j,ii,jj;
    for(i=0; i<=logn; i++)
        for(j=0; j<=logm; j++) {
            if(i==0&&j==0) continue;
            for(ii=1; ii+p[i]<=n+1; ii++)
                for(jj=1; jj+p[j]<=m+1; jj++) {
                    if(i==0)
                        dp[i][j][ii][jj]=max(dp[i][j-1][ii][jj],dp[i][j-1][ii][jj+p[j-1]]);
                    else
                        dp[i][j][ii][jj]=max(dp[i-1][j][ii][jj],dp[i-1][j][ii+p[i-1]][jj]);
                }
        }
}
bool check() {
    if(ans==dp[0][0][r1][c1]) return 1;
    if(ans==dp[0][0][r1][c2]) return 1;
    if(ans==dp[0][0][r2][c1]) return 1;
    if(ans==dp[0][0][r2][c2]) return 1;
    return 0;
}
int main() {
    int i,j,q;
    while(~scanf("%d%d",&n,&m)) {
        for(i=1; i<=n; i++)
            for(j=1; j<=m; j++)
                scanf("%d",&dp[0][0][i][j]);
        redeal();
        scanf("%d",&q);
        while(q--) {
            scanf("%d%d%d%d",&r1,&c1,&r2,&c2);
            int log1=mylog(r2-r1+1);
            int log2=mylog(c2-c1+1);
            int t1=dp[log1][log2][r1][c1];
            int t2=dp[log1][log2][r2-p[log1]+1][c1];
            int t3=dp[log1][log2][r1][c2-p[log2]+1];
            int t4=dp[log1][log2][r2-p[log1]+1][c2-p[log2]+1];
            ans=max(max(t1,t2),max(t3,t4));
            printf("%d ",ans);
            if(check()) printf("yes\n");
            else printf("no\n");
        }
    }
}
//划分树（第K小数）
int val[100010][20],L[100010][20],sorted[100010];
void creat(int pos,int cur,int l,int r) {
    int i;
    if(l==r) return ;
    int lpos,rpos,mid=(l+r)>>1,lsame,same;
    lsame=mid+1-l;
    for(i=l; i<=r; i++) if(val[i][cur]<sorted[mid]) lsame--;
    same=0,lpos=l,rpos=mid+1;
    for(i=l; i<=r; i++) {
        if(i==l) L[i][cur]=0;
        else L[i][cur]=L[i-1][cur];
        if(val[i][cur]<sorted[mid]) val[lpos++][cur+1]=val[i][cur],L[i][cur]++;
        else if(val[i][cur]>sorted[mid]) val[rpos++][cur+1]=val[i][cur];
        else {
            if(lsame>same) {
                same++;
                L[i][cur]++;
                val[lpos++][cur+1]=val[i][cur];
            } else val[rpos++][cur+1]=val[i][cur];
        }
    }
    creat(pos*2,cur+1,l,mid);
    creat(pos*2+1,cur+1,mid+1,r);
}
int sel(int pos,int cur,int l,int r,int start,int end,int k) {
    if(start==end) return val[start][cur];
    int s,ss;
    if(start==l) {
        s=L[end][cur];
        ss=0;
    } else {
        s=L[end][cur]-L[start-1][cur];
        ss=L[start-1][cur];
    }
    int mid=(l+r)>>1;
    if(s>=k) return sel(pos*2,cur+1,l,mid,l+ss,l+ss+s-1,k);
    else {
        int sb=start-l-ss;
        int sa=end+1-start-s;
        return sel(pos*2+1,cur+1,mid+1,r,mid+1+sb,mid+sa+sb,k-s);
    }
}
int main() {
    int n,m,i,a,b,c;
    scanf("%d%d",&n,&m);
    for(i=1; i<=n; i++) {
        scanf("%d",&sorted[i]);
        val[i][0]=sorted[i];
    }
    sort(sorted+1,sorted+n+1);
    creat(1,0,1,n);
    while(m--) {
        scanf("%d%d%d",&a,&b,&c);
        printf("%d\n",sel(1,0,1,n,a,b,c));
    }
    return 0;
}
//字符串循环最小
int main() {
    int T,t,i,j,k,st,len;
    scanf("%d",&T);
    while(T--) {
        scanf("%s",str);
        len=strlen(str);
        i=0,j=1,k=0;
        while(i<len&&j<len&&k<len) {
            t=str[(i+k)%len]-str[(j+k)%len];
            if(!t) k++;
            else {
                if(t>0) i+=k+1;
                if(t<0) j+=k+1;
                if(i==j) j++;
                k=0;
            }
        }
        st=min(i,j);
        for(i=st; i<st+len; i++) printf("%c",str[i%len]);
        printf("\n");
    }
    return 0;
}
//LAZY操作线段树zju 3379
const double pi=acos(-1.0),eps=1e-6;
int sgn(double x) {
    return (x>eps)-(x<-eps);
}
double cal(double x,double y,double a) {
    double len=sqrt(x*x+y*y);
    double s=(sqrt(1+4*a*a*(x*x+y*y))-1.0)/(2*a*a);
    double tx=sqrt(s);
    return asin(tx/len);
}
struct node {
    double s,t;
    node() {}
    node(double a,double b):s(a),t(b) {}
} rec[60010];
double A[120010];
struct gg {
    int w,lazy;
} T[500010];
void creat(int pos,int l,int r) {
    T[pos].w=T[pos].lazy=0;
    if(l==r) return ;
    int mid=(l+r)>>1;
    creat(pos*2,l,mid);
    creat(pos*2+1,mid+1,r);
}
void update(int pos) {
    T[pos].w=max(T[pos*2].w,T[pos*2+1].w);
}
void insert(int pos,int l,int r,int start,int end) {
    if(l==start&&r==end) {
        T[pos].w++;
        T[pos].lazy++;
        return ;
    }
    int mid=(l+r)>>1;
    int tmp=T[pos].lazy;
    T[pos*2].lazy+=tmp;
    T[pos*2].w+=tmp;
    T[pos*2+1].lazy+=tmp;
    T[pos*2+1].w+=tmp;
    T[pos].lazy=0;
    if(end<=mid) insert(pos*2,l,mid,start,end);
    else if(start>mid) insert(pos*2+1,mid+1,r,start,end);
    else {
        insert(pos*2,l,mid,start,mid);
        insert(pos*2+1,mid+1,r,mid+1,end);
    }
    update(pos);
}
int m;
int find(double a) {
    int low=0,high=m-1,mid;
    while(low<=high) {
        mid=(low+high)>>1;
        if(sgn(a-A[mid])==0) return mid;
        if(a>A[mid]) low=mid+1;
        else high=mid-1;
    }
}
double xx[30010],yy[30010];
int main() {
    int n,i,cnt,ss,ee;
    double a,start,thera,x,y,s,t;
    while(~scanf("%d%lf",&n,&a)) {
        if(!n) {
            printf("0 daze\n");
            continue;
        }
        cnt=0;
        for(i=0; i<n; i++) scanf("%lf",&xx[i]);
        for(i=0; i<n; i++) scanf("%lf",&yy[i]);
        for(i=0; i<n; i++) {
            x=xx[i],y=yy[i];
            thera=cal(x,y,a);
            start=acos(x/(sqrt(x*x+y*y)));
            if(sgn(y)<0) start=2*pi-start;
            s=start-thera,t=start+thera;
            if(sgn(s)<0) rec[cnt++]=node(2*pi+s,pi*2),rec[cnt++]=node(0,t);
            else if(sgn(t-2*pi)>0) rec[cnt++]=node(s,2*pi),rec[cnt++]=node(0,t-2*pi);
            else {
                rec[cnt++]=node(s,t);
                if(sgn(s)==0) rec[cnt++]=node(pi*2,pi*2);
                if(sgn(t-2*pi)==0) rec[cnt++]=node(0,0);
            }
        }
        m=cnt<<1;
        for(i=0; i<cnt; i++) A[i<<1]=rec[i].s,A[(i<<1)+1]=rec[i].t;
        sort(A,A+m);
        m=unique(A,A+m)-A;
        creat(1,0,m-1);
        for(i=0; i<cnt; i++) {
            ss=find(rec[i].s);
            ee=find(rec[i].t);
            insert(1,0,m-1,ss,ee);
        }
        printf("%d daze\n",T[1].w);
    }
    return 0;
}
//多维几何体求交
const ll M=14121413LL;
struct node {
    int bot[10],top[10];
} C[100010],D[100010];
int n,k,cnt;
bool inter(node a,node b) {
    int i;
    for(i=0; i<n; i++) if(a.bot[i]>=b.top[i]||a.top[i]<=b.bot[i]) return 0;
    return 1;
}
void ins(node a) {
    int i;
    for(i=0; i<n; i++) if(a.bot[i]>=a.top[i]) return ;
    D[k++]=a;
}
void broke(node a,node b) {
    node temp=a,tt;
    int i;
    for(i=0; i<n; i++) {
        if(b.bot[i]>=temp.bot[i]&&b.bot[i]<=temp.top[i]) {
            tt=temp;
            tt.top[i]=b.bot[i];
            ins(tt);
            temp.bot[i]=b.bot[i];
        }
        if(b.top[i]>=temp.bot[i]&&b.top[i]<=temp.top[i]) {
            tt=temp;
            tt.bot[i]=b.top[i];
            ins(tt);
            temp.top[i]=b.top[i];
        }
    }
}
void insert(node c) {
    int i;
    k=0;
    for(i=0; i<cnt; i++)
        if(inter(C[i],c)) broke(C[i],c);
        else D[k++]=C[i];
    D[k++]=c;
    for(i=0; i<k; i++) C[i]=D[i];
    cnt=k;
}
void sum() {
    ll k,ans=0;
    int i,j;
    for(i=0; i<cnt; i++) {
        k=1;
        for(j=0; j<n; j++) k=(k*(ll)(C[i].top[j]-C[i].bot[j]))%M;
        ans=(ans+k)%M;
    }
    printf("%I64d\n",ans);
}
int main() {
    node temp;
    int m,i;
    while(~scanf("%d%d",&m,&n)) {
        cnt=0;
        while(m--) {
            for(i=0; i<n; i++) scanf("%d",&temp.bot[i]);
            for(i=0; i<n; i++) scanf("%d",&temp.top[i]);
            insert(temp);
        }
        sum();
    }
    return 0;
}
//DLX重复覆盖 hdu 2295
const int inf=0x3fffffff,MAXN=3000;
const double eps=1e-8;
int op,m,ans,n,s,l;
bool hash[70];
struct node {
    double x,y;
} city[60],radar[60];
int L[MAXN],R[MAXN],D[MAXN],U[MAXN],S[MAXN],C[MAXN];
void init() {
    int i;
    for(i=0; i<=m; i++)
        L[i]=i-1,R[i]=i+1,U[i]=i,D[i]=i,S[i]=0;
    L[0]=m,R[m]=0;
    op=m+1;
}
void insert(int head,int &pre,int pos) {
    U[op]=U[pos],D[U[pos]]=op;
    D[op]=pos,U[pos]=op;
    S[pos]++,C[op]=pos;
    L[op]=pre,R[pre]=op;
    R[op]=head,L[head]=op;
    pre=op++;
}
void remove(int c) {
    int i;
    for(i=D[c]; i!=c; i=D[i]) L[R[i]]=L[i],R[L[i]]=R[i];
}
void resume(int c) {
    int i;
    for(i=U[c]; i!=c; i=U[i]) L[R[i]]=i,R[L[i]]=i;
}
int h() {
    int ret=0,i,j,l;
    for(i=0; i<m; i++) hash[i]=0;
    for(i=R[m]; i!=m; i=R[i]) {
        if(!hash[i]) {
            ret++;
            for(j=D[i]; j!=i; j=D[j])
                for(l=R[j]; l!=j; l=R[l]) hash[C[l]]=1;
        }
    }
    return ret;
}
bool dlx(int k,int lim) {
    if(k+h()>lim) return 0;
    if(R[m]==m) {
        return 1;
    }
    int c,mins=inf,i,j;
    for(i=R[m]; i!=m; i=R[i]) if(mins>S[i]) c=i,mins=S[i];
    for(i=D[c]; i!=c; i=D[i]) {
        remove(i);
        for(j=R[i]; j!=i; j=R[j]) remove(j);
        if(dlx(k+1,lim)) return 1;
        for(j=L[i]; j!=i; j=L[j]) resume(j);
        resume(i);
    }
    return 0;
}
double sqr(double x) {
    return x*x;
}
double dis(node p1,node p2) {
    return sqrt(sqr(p1.x-p2.x)+sqr(p1.y-p2.y));
}
int sgn(double x) {
    return (x>eps)-(x<-eps);
}
void addedge(double r) {
    int i,j,head,pre;
    for(i=0; i<l; i++) {
        head=pre=op;
        for(j=0; j<n; j++) if(sgn(dis(city[j],radar[i])-r)<=0) insert(head,pre,j);
    }
}
void solve() {
    double low=0,high=10000,mid;
    while(high-low>eps) {
        mid=(low+high)/2;
        init();
        addedge(mid);
        if(dlx(0,s)) high=mid;
        else low=mid;
    }
    printf("%.6lf\n",high);
}
int main() {
    int T,i;
    scanf("%d",&T);
    while(T--) {
        scanf("%d%d%d",&n,&l,&s);
        m=n;
        for(i=0; i<n; i++) scanf("%lf%lf",&city[i].x,&city[i].y);
        for(i=0; i<l; i++) scanf("%lf%lf",&radar[i].x,&radar[i].y);
        solve();
    }
    return 0;
}
//DLX精确覆盖
int L[MAXN],R[MAXN],D[MAXN],U[MAXN],S[MAXN],C[MAXN];
void init() {
    int i;
    for(i=0; i<=m; i++)
        L[i]=i-1,R[i]=i+1,U[i]=i,D[i]=i,S[i]=0;
    L[0]=m,R[m]=0;
    op=m+1;
}
void insert(int head,int &pre,int pos) {
    U[op]=U[pos],D[U[pos]]=op;
    D[op]=pos,U[pos]=op;
    S[pos]++,C[op]=pos;
    L[op]=pre,R[pre]=op;
    R[op]=head,L[head]=op;
    pre=op++;
}
void remove(int c) {
    L[R[c]]=L[c];
    R[L[c]]=R[c];
    int i,j;
    for(i=D[c]; i!=c; i=D[i])
        for(j=R[i]; j!=i; j=R[j]) {
            U[D[j]]=U[j];
            D[U[j]]=D[j];
            S[C[j]]--;
        }
}
void resume(int c) {
    int i,j;
    for(i=U[c]; i!=c; i=U[i])
        for(j=L[i]; j!=i; j=L[j]) {
            U[D[j]]=j;
            D[U[j]]=j;
            S[C[j]]++;
        }
    L[R[c]]=c;
    R[L[c]]=c;
}
void dlx(int k) {
    if(k>=ans) return ;
    if(R[m]==m) {
        ans=k;
        return;
    }
    int c,mins=inf,i,j;
    for(i=R[m]; i!=m; i=R[i]) if(mins>S[i]) c=i,mins=S[i];
    remove(c);
    for(i=D[c]; i!=c; i=D[i]) {
        for(j=R[i]; j!=i; j=R[j]) remove(C[j]);
        dlx(k+1);
        for(j=L[i]; j!=i; j=L[j]) resume(C[j]);
    }
    resume(c);
}
void solve() {
    int i;
    for(i=0; i<m; i++) if(!S[i]) {
            printf("-1\n");
            return;
        }
    ans=inf;
    dlx(0);
    if(ans!=inf) printf("%d\n",ans);
    else printf("-1\n");
}
int main() {
    int T,r,c,p,head,pre,x1,x2,y1,y2,i,j;
    scanf("%d",&T);
    while(T--) {
        scanf("%d%d%d",&r,&c,&p);
        m=r*c;
        init();
        while(p--) {
            head=op,pre=op;
            scanf("%d%d%d%d",&x1,&y1,&x2,&y2);
            for(i=x1; i<x2; i++)
                for(j=y1; j<y2; j++)
                    insert(head,pre,i*c+j);
        }
        solve();
    }
    return 0;
}
//带权并查集
struct node {
    int father,rel;
} s[50010];
int find(int a) {
    if(s[a].father==a) return a;
    int t=find(s[a].father);
    s[a].rel=(s[a].rel+s[s[a].father].rel)%3;
    return s[a].father=t;
}
int main() {
    int n,k,i,ans=0,d,a,b,ta,tb;
    scanf("%d%d",&n,&k);
    for(i=1; i<=n; i++) s[i].father=i,s[i].rel=0;
    while(k--) {
        scanf("%d%d%d",&d,&a,&b);
        if(a>n||b>n) {
            ans++;
            continue;
        }
        d--;
        ta=find(a),tb=find(b);
        if(ta==tb) {
            if((s[a].rel-s[b].rel+3)%3!=d) {
                ans++;
                continue;
            }
        }
        s[ta].father=tb;
        s[ta].rel=(s[b].rel-s[a].rel+d+3)%3;
    }
    printf("%d\n",ans);
    return 0;
}
//Splay
const int Size=200010,Inf=0x3fffffff;
struct node {
    int v,pre,s,rev,m,add;
    int ch[2];
} T[Size];
int n,root,op;
void init() {
    op=1;
    int i;
    T[op].v=Inf,T[op].pre=op+1,T[op].ch[0]=T[op].ch[1]=-1;
    T[op].s=1,T[op].rev=0,T[op].m=Inf,T[op].add=0;
    for(i=1; i<=n+1; i++) {
        op++;
        if(i!=n+1) scanf("%d",&T[op].v);
        else T[op].v=Inf;
        T[op].pre=op+1,T[op].ch[0]=op-1,T[op].ch[1]=-1;
        T[op].s=i+1,T[op].rev=0,T[op].m=min(T[op-1].m,T[op].v),T[op].add=0;
    }
    T[op].pre=0;
    root=op;
}
int sz(int x) {
    if(~x) return T[x].s;
    return 0;
}
void update(int x) {
    int &l=T[x].ch[0],&r=T[x].ch[1];
    if(T[x].rev) {
        T[x].rev=0;
        if(~l) T[l].rev^=1;
        if(~r) T[r].rev^=1;
        swap(l,r);
    }
    if(~l) {
        T[l].add+=T[x].add;
        T[l].m+=T[x].add;
        T[l].v+=T[x].add;
    }
    if(~r) {
        T[r].add+=T[x].add;
        T[r].m+=T[x].add;
        T[r].v+=T[x].add;
    }
    T[x].add=0;
}
int rank(int k,int x) {
    update(x);
    int ln=sz(T[x].ch[0]);
    if(ln+1==k) return x;
    if(k<=ln) return rank(k,T[x].ch[0]);
    if(k>ln+1) return rank(k-ln-1,T[x].ch[1]);
}
int find(int x) {
    if(~x) return T[x].m;
    return Inf;
}
void tran(int x) {
    T[x].s=sz(T[x].ch[0])+sz(T[x].ch[1])+1;
    T[x].m=min(T[x].v,min(find(T[x].ch[0]),find(T[x].ch[1])));
}
void R(int x,int judge) {
    int f=T[x].pre;
    int gf=T[f].pre;
    T[gf].ch[T[gf].ch[1]==f]=x;
    T[f].pre=x;
    T[f].ch[!judge]=T[x].ch[judge];
    tran(f);
    T[T[x].ch[judge]].pre=f;
    T[x].pre=gf;
    T[x].ch[judge]=f;
    tran(x);
}
void splay(int x,int r) {
    int f,gf;
    while(T[x].pre!=r) {
        f=T[x].pre;
        if(T[f].pre==r) R(x,T[f].ch[0]==x);
        else {
            gf=T[f].pre;
            if((T[gf].ch[0]==f)^(T[f].ch[0]==x)) R(x,T[f].ch[0]==x),R(x,T[gf].ch[0]==x);
            else R(f,T[gf].ch[0]==f),R(x,T[f].ch[0]==x);
        }
    }
    if(r==0) root=x;
}
void up(int a) {
    int x=a;
    while(x) {
        tran(x);
        x=T[x].pre;
    }
}
void Reverse(int l,int r) {
    int a=rank(l-1,root),b=rank(r+1,root);
    splay(b,0);
    splay(a,b);
    T[T[a].ch[1]].rev^=1;
}
void Add(int l,int r,int w) {
    int a=rank(l-1,root),b=rank(r+1,root);
    splay(b,0);
    splay(a,b);
    T[T[a].ch[1]].add+=w;
    T[T[a].ch[1]].v+=w;
    T[T[a].ch[1]].m+=w;
}
int Min(int l,int r) {
    int a=rank(l-1,root),b=rank(r+1,root);
    splay(b,0);
    splay(a,b);
    return T[T[a].ch[1]].m;
}
void Insert(int x,int w) {
    int a=rank(x,root),b=rank(x+1,root);
    splay(b,0);
    splay(a,b);
    op++;
    T[op].v=w,T[op].pre=b,T[op].s=T[a].s+1,T[op].rev=0;
    T[op].m=min(T[a].m,w),T[op].add=0,T[op].ch[0]=a,T[op].ch[1]=-1;
    T[a].pre=op;
    T[b].ch[0]=op;
    tran(b);
}
void Delete(int x) {
    int a=rank(x,root),b=rank(x+1,root);
    splay(b,0);
    splay(a,b);
    T[b].ch[0]=T[a].ch[0];
    tran(b);
    if(~T[b].ch[0]) T[T[b].ch[0]].pre=b;
}
int Left(int x) {
    update(x);
    if(T[x].ch[0]==-1) return x;
    return Left(T[x].ch[0]);
}
void Revolve(int l,int r,int w) {
    int c=(w%(r+1-l)+(r+1-l))%(r+1-l);
    if(c==0||r==l) return ;
    int a=l,b=r-c;
    c=r;
    a=rank(a-1,root);
    b=rank(b,root);
    c=rank(c+1,root);
    splay(b,0);
    splay(a,b);
    splay(c,b);
    int d=T[a].ch[1],e=T[c].ch[0];
    T[c].ch[0]=-1;
    tran(c);
    if(~d) {
        int f=Left(d);
        T[a].ch[1]=e,T[a].pre=f;
        T[e].pre=a;
        T[f].ch[0]=a;
        T[d].pre=b;
        T[b].ch[0]=d;
    } else {
        T[a].ch[1]=e;
        T[e].pre=a;
    }
    up(a);
}
int main() {
    int m,l,r,w;
    char str[10];
    while(~scanf("%d",&n)) {
        init();
        scanf("%d",&m);
        while(m--) {
            scanf("%s",str);
            if(!strcmp(str,"ADD")) {
                scanf("%d%d%d",&l,&r,&w);
                Add(l+1,r+1,w);
            }
            if(!strcmp(str,"REVERSE")) {
                scanf("%d%d",&l,&r);
                Reverse(l+1,r+1);
            }
            if(!strcmp(str,"INSERT")) {
                scanf("%d%d",&l,&w);
                Insert(l+1,w);
            }
            if(!strcmp(str,"DELETE")) {
                scanf("%d",&l);
                Delete(l+1);
            }
            if(!strcmp(str,"REVOLVE")) {
                scanf("%d%d%d",&l,&r,&w);
                Revolve(l+1,r+1,w);
            }
            if(!strcmp(str,"MIN")) {
                scanf("%d%d",&l,&r);
                printf("%d\n",Min(l+1,r+1));
            }
        }
    }
    return 0;
}
//动态树 Hdu 2475
#include<iostream>
#include<stdio.h>
using namespace std;
const int Size=50010;
struct node {
    int pre;
    int ch[2];
} T[Size];
void R(int x,int judge) {
    int f=T[x].pre;
    int gf=T[f].pre;
    if(gf!=-1) {
        if(T[gf].ch[0]==f) T[gf].ch[0]=x;
        if(T[gf].ch[1]==f) T[gf].ch[1]=x;
    }
    T[f].pre=x;
    T[f].ch[!judge]=T[x].ch[judge];

    T[T[x].ch[judge]].pre=f;
    T[x].pre=gf;
    T[x].ch[judge]=f;
}
void splay(int x,int r) {
    int f,gf;
    while(T[x].pre!=r) {
        f=T[x].pre;
        if(T[f].pre==r) R(x,T[f].ch[0]==x);
        else {
            gf=T[f].pre;
            if((T[gf].ch[0]==f)^(T[f].ch[0]==x)) R(x,T[f].ch[0]==x),R(x,T[gf].ch[0]==x);
            else R(f,T[gf].ch[0]==f),R(x,T[f].ch[0]==x);
        }
    }
}
int find(int x) {
    int tmp=T[x].pre,tt=x;
    while((~tmp)&&(T[tmp].ch[0]==tt||T[tmp].ch[1]==tt)) tt=tmp,tmp=T[tmp].pre;
    return tmp;
}
void splice(int x) {
    int p=T[x].pre;
    T[p].ch[0]=x;
}
void expose(int x) {
    int p,tmp;
    while(1) {
        tmp=find(x);
        splay(x,tmp);
        p=T[x].pre;
        if(p==-1) break;
        tmp=find(p);
        splay(p,tmp);
        splice(x);
    }
    T[x].ch[0]=-1;
}
int root(int x) {
    expose(x);
    int t=x;
    while(T[t].ch[1]!=-1) t=T[t].ch[1];
    return t;
}
int parent(int x) {
    expose(x);
    int t=T[x].ch[1];
    if(t==-1) return -1;
    while(T[t].ch[0]!=-1) t=T[t].ch[0];
    return t;
}
void link(int a,int b) {
    expose(a);
    T[a].pre=b;
}
void cut(int a) {
    expose(a);
    int p=T[a].ch[1];
    if(p!=-1) T[p].pre=-1;
    T[a].ch[1]=-1;
}
int main() {
    int i,a,b,q,n,cases=0;
    char str[20];
    while(~scanf("%d",&n)) {
        for(i=0; i<n; i++) T[i].pre=T[i].ch[0]=T[i].ch[1]=-1;
        for(i=0; i<n; i++) {
            scanf("%d",&a);
            if(a) link(i,a-1);
        }
        if(cases) puts("");
        cases++;
        scanf("%d",&q);
        while(q--) {
            scanf("%s",str);
            if(str[0]=='Q') {
                scanf("%d",&a);
                a--;
                printf("%d\n",root(a)+1);
            } else {
                scanf("%d%d",&a,&b);
                a--,b--;
                int c=parent(a);
                cut(a);
                if(~b) {
                    if(root(b)==a) {
                        if(~c) link(a,c);
                    } else link(a,b);
                }
            }
        }
    }
    return 0;
}
