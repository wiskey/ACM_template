//Miller-rabin素数测试
__int64 mod(__int64 a,__int64 m,__int64 p) {
    __int64 ans=1,power=a%p,temp=m;
    while(temp) {
        if(temp&1) ans=(ans*power)%p;
        power=(power*power)%p;
        temp>>=1;
    }
    return ans;
}
bool check(__int64 a,__int64 n) {
    __int64 len=0,m=n-1,i;
    while(m&&(!(m&1))) {
        m>>=1;
        len++;
    }
    __int64 b=mod(a,m,n);
    if(b==1||b==n-1) return true;
    for(i=1; i<len; i++) {
        b=(b*b)%n;
        if(b==n-1) return true;
    }
    return false;
}
bool prime(__int64 n) {
    __int64 i,a;
    if(n==1) return false;
    if(n==2) return true;
    for(i=1; i<=5; i++) {
        a=rand()%(n-2)+2;
        if(!check(a,n))
            return false;
    }
    return true;
}
int main() {
    __int64 n;
    srand(time(NULL));
    while(~scanf("%I64d",&n)) printf("%d\n",prime(n));

}
线性筛素数表
int ans[3000010]= {0};
bool mark[3000010]= {0};
int rec[3000010];
void init() {
    int size=0,i,j;
    mark[1]=1;
    for(i=2; i<=n; i++) {
        if(!mark[i]) rec[size++]=i;
        for(j=0; (j<size&&i*rec[j]<=n); j++) mark[i*rec[j]]=1;
    }
}
扩展欧几里得
int extgcd(int a,int b,int &x,int &y) {
    if(b==0) {
        x=1;
        y=0;
        return a;
    }
    int d=extgcd(b,a%b,x,y);
    int temp=x;
    x=y;
    y=temp-a/b*y;
    return d;
}
int main() {
    int lim,a,b,c,k,x,y,d;
    while(scanf("%d%d%d%d",&a,&b,&c,&k)==4) {
        lim=1<<k;
        d=extgcd(lim,c,x,y);
        b-=a;
        if(b%d) printf("FOREVER\n");
        else {
            y*=(b/d);
            lim/=d;
            printf("%d\n",(y%lim+lim)%lim);
        }
    }
}
高斯消元
double rec[110][110];
double x[110];
bool sure_x[110];
int n,m;
int solve() {
    int cur,i,j,row,col,l,num;
    double temp;
    for(row=0,col=0; row<n&&col<m; row++,col++) {
        cur=row;
        for(j=row+1; j<n; j++) if(fabs(rec[j][col])>fabs(rec[cur][col])) cur=j;
        if(rec[cur][col]==0) {
            row--;
            continue;
        }
        for(j=col; j<=m; j++) swap(rec[row][j],rec[cur][j]);
        for(j=row+1; j<n; j++) {
            temp=rec[j][col]/rec[row][col];
            for(l=col; l<=m; l++) rec[j][l]-=rec[row][l]*temp;
        }
    }
    memset(sure_x,0,sizeof(sure_x));
    for(i=row; i<n; i++) if(rec[i][m]!=0) return -1;
    if(row<m) {
        for(i=row-1; i>=0; i--) {
            num=0;
            for(j=0; j<m; j++) if(rec[i][j]!=0&&!sure_x[j]) num++,cur=j;
            if(num>1) continue;
            temp=rec[i][m];
            for(j=0; j<m; j++) if(rec[i][j]!=0&&j!=cur) temp-=rec[i][j]*x[j];
            x[cur]=temp/rec[i][cur];
            sure_x[cur]=1;
        }
        return m-row;
    }
    for(i=n-1; i>=0; i--) {
        temp=rec[i][m];
        for(j=i+1; j<m; j++) if(rec[i][j]!=0) temp-=rec[i][j]*x[j];
        x[i]=temp/rec[i][i];
    }
    return 0;
}
int main() {
    int i,j;
    scanf("%d%d",&n,&m);
    for(i=0; i<n; i++)
        for(j=0; j<=m; j++)
            scanf("%lf",&rec[i][j]);
    int judge=solve();
    if(judge==-1) printf("no solution");
    if(judge==0) for(i=0; i<m; i++) printf("%.2lf ",x[i]);
    if(judge!=0&&judge!=-1) for(i=0; i<m; i++) {
            if(sure_x[i]) printf("%.2lf ",x[i]);
            else printf("free ");
        }
    printf("\n");
}
辗转相除的高斯消元
ll rec[210][210];
int n;
ll p;
void cut(ll &x,ll &y,int ss,ll *a,ll *b) {
    int i;
    if(y==0) return ;
    for(i=ss; i<n; i++) a[i]-=x/y*b[i],a[i]=(a[i]%p+p)%p;
    x%=y;
    cut(y,x,ss,b,a);
}
void gauss() {
    int row,col,i,l,j,pp=0;
    ll ans;
    for(row=0,col=0; row<n&&col<n; row++,col++) {
        for(i=row; i<n; i++) if(rec[i][col]) break;
        if(i==n) {
            printf("0\n");
            return ;
        }
        for(j=i+1; j<n; j++) {
            if(rec[j][col]==0) continue;
            cut(rec[i][col],rec[j][col],col+1,rec[i],rec[j]);
            if(rec[i][col]==0) {
                pp++;
                for(l=col; l<n; l++) swap(rec[i][l],rec[j][l]);
            }
        }
        if(row!=i) pp++;
        for(j=col; j<n; j++) swap(rec[row][j],rec[i][j]);
    }
    if(pp&1) ans=-1;
    else ans=1;
    for(i=0; i<n; i++) ans*=rec[i][i],ans=(ans%p+p)%p;
    printf("%I64d\n",ans);
}
int main() {
    int i,j;
    while(~scanf("%d%I64d",&n,&p)) {
        for(i=0; i<n; i++)
            for(j=0; j<n; j++)
                scanf("%I64d",&rec[i][j]);
        gauss();
    }
    return 0;
}
Polya+欧拉函数+矩阵法
bool mark[Max+1];
int rec[Max],num[Max],seq[Max];
int matrix[11][11];
int size,n,m,res,cnt;
void init() {
    int i,j;
    size=0;
    for(i=1; i<=Max; i++) mark[i]=1;
    for(i=2; i<=Max; i++) {
        if(mark[i]) rec[size++]=i;
        for(j=0; j<size&&rec[j]*i<=Max; j++) mark[rec[j]*i]=0;
    }
}
int phi(int num) {
    int tn=num,i=0,ans=num;
    while(1) {
        if(rec[i]*rec[i]>tn) break;
        if(tn%rec[i]==0) {
            ans-=ans/rec[i];
            while(tn%rec[i]==0) tn/=rec[i];
        }
        i++;
    }
    if(tn!=1) ans-=ans/tn;
    return ans;
}
void mult(int at[11][11],int bt[11][11]) {
    int temp[11][11]= {0};
    int i,j,k;
    for(i=1; i<=m; i++)
        for(j=1; j<=m; j++) {
            for(k=1; k<=m; k++)
                temp[i][j]=(temp[i][j]+at[i][k]*bt[k][j]);
            temp[i][j]%=mod;
        }
    for(i=1; i<=m; i++)
        for(j=1; j<=m; j++) at[i][j]=temp[i][j];
}
int deal(int num) {
    int sum,temp=num;
    int ans[11][11],power[11][11];
    int i,j;
    for(i=1; i<=m; i++)
        for(j=1; j<=m; j++) {
            if(i==j) ans[i][j]=1;
            else ans[i][j]=0;
            power[i][j]=matrix[i][j];
        }
    while(temp) {
        if(temp&1) mult(ans,power);
        mult(power,power);
        temp>>=1;
    }
    sum=0;
    for(i=1; i<=m; i++) sum=(sum+ans[i][i])%mod;
    return sum;
}
void dfs(int pos,int w) {
    if(pos==cnt) {
        res=(res+phi(n/w)%mod*(deal(w)%mod))%mod;
        return ;
    }
    int sum=1,i;
    for(i=0; i<=num[pos]; i++) {
        dfs(pos+1,w*sum);
        sum*=seq[pos];
    }
}
void gcd(int a,int b,int &x,int &y) {
    if(b==0) {
        x=1,y=0;
        return ;
    }
    gcd(b,a%b,x,y);
    int tmp=x;
    x=y;
    y=tmp-a/b*x;
}
int get_inverse() {
    int x,y;
    gcd(n,mod,x,y);
    return (x%mod+mod)%mod;
}
void solve() {
    int i,tn=n,sum=1;
    cnt=0;
    for(i=0; i<=size; i++) num[i]=0;
    i=0;
    while(1) {
        if(rec[i]*rec[i]>tn) break;
        while(tn%rec[i]==0) tn/=rec[i],num[cnt]++,seq[cnt]=rec[i];
        if(num[cnt]) cnt++;
        i++;
    }
    if(tn!=1) num[cnt]++,seq[cnt]=tn,cnt++;
    res=0;
    for(i=0; i<=num[0]; i++) {
        dfs(1,sum);
        sum*=seq[0];
    }
    int inv=get_inverse();
    printf("%d\n",res*inv%mod);
}
int main() {
    int T,k,i,a,b,j;
    init();
    scanf("%d",&T);
    while(T--) {
        scanf("%d%d%d",&n,&m,&k);
        for(i=1; i<=m; i++) for(j=1; j<=m; j++) matrix[i][j]=1;
        while(k--) {
            scanf("%d%d",&a,&b);
            matrix[a][b]=matrix[b][a]=0;
        }
        if(n==1) {
            printf("%d\n",m%mod);
            continue;
        }
        solve();
    }
    return 0;
}
Java分数高斯消元
import java.io.*;
import java.util.*;
import java.math.*;
public class Main {
    class node {
        boolean judge;
        BigInteger up,down;
        node(boolean a,BigInteger b,BigInteger c) {
            judge=a;
            up=b;
            down=c;
        }
        int comp(node a) {
            return up.multiply(a.down).compareTo(a.up.multiply(down));
        }
        node tim(node a) {
            if(a.up.equals(zero)||up.equals(zero)) return new node(true,zero,BigInteger.valueOf(1));
            BigInteger n1=up.multiply(a.up);
            BigInteger n2=down.multiply(a.down);
            BigInteger n3=n1.gcd(n2);
            return new node(a.judge==judge,n1.divide(n3),n2.divide(n3));
        }
        node min(node a) {
            if(up.equals(zero)) return new node(!a.judge,a.up,a.down);
            if(a.up.equals(zero)) return new node(judge,up,down);
            if((!judge&&a.judge)||(judge&&!a.judge)) {
                BigInteger n1=up.multiply(a.down).add(a.up.multiply(down));
                BigInteger n2=down.multiply(a.down);
                BigInteger n3=n1.gcd(n2);
                return new node(judge,n1.divide(n3),n2.divide(n3));
            }
            BigInteger n1=up.multiply(a.down).subtract(a.up.multiply(down));
            BigInteger n2=down.multiply(a.down);
            if(n1.equals(zero)) return new node(true,zero,BigInteger.valueOf(1));
            boolean p=n1.compareTo(zero)>0;
            BigInteger n3=n1.abs().gcd(n2);
            return new node(judge==p,n1.abs().divide(n3),n2.divide(n3));
        }
        node div(node a) {
            if(up.equals(zero)) return new node(true,zero,BigInteger.valueOf(1));
            BigInteger n1=up.multiply(a.down);
            BigInteger n2=down.multiply(a.up);
            BigInteger n3=n1.gcd(n2);
            return new node(a.judge==judge,n1.divide(n3),n2.divide(n3));
        }
        void print() {
            if(up.equals(zero)) {
                System.out.println(0);
                return ;
            }
            if(!judge) System.out.print("-");
            System.out.print(up);
            if(!down.equals(BigInteger.valueOf(1))) System.out.print("/"+down);
            System.out.println();
        }
    };
    int n;
    BigInteger zero=BigInteger.valueOf(0);
    node[][] a=new node[110][110];
    void gauss() {
        int i,j,k,l,p;
        node[] x=new node[110];
        for(i=0,k=0; i<n&&k<n; i++,k++) {
            p=i;
            for(j=i+1; j<n; j++) if(a[j][k].comp(a[p][k])>0) p=j;
            if(a[p][k].up.equals(zero)) {
                i--;
                continue;
            }
            for(j=k; j<=n; j++) {
                node b=a[i][j];
                a[i][j]=a[p][j];
                a[p][j]=b;
            }
            for(j=i+1; j<n; j++) {

                if(a[j][k].up.equals(zero)) continue;
                node b=a[j][k].div(a[i][k]);
                for(l=k; l<=n; l++) a[j][l]=a[j][l].min(b.tim(a[i][l]));
            }
        }
        k=i;
        for(i=k; i<n; i++) if(!a[i][n].up.equals(zero)) {
                System.out.println("No solution.");
                System.out.println();
                return ;
            }
        for(i=n-1; i>=0; i--) {
            node temp=a[i][n];
            for(j=i+1; j<n; j++) temp=temp.min(a[i][j].tim(x[j]));
            x[i]=temp.div(a[i][i]);
        }
        for(i=0; i<n; i++) x[i].print();
        System.out.println();
    }
    public Main() {
        Scanner in=new Scanner(System.in);
        int i,j;
        BigInteger temp;
        while(in.hasNext()) {
            n=in.nextInt();
            for(i=0; i<n; i++)
                for(j=0; j<=n; j++) {
                    temp=in.nextBigInteger();
                    a[i][j]=new node(temp.compareTo(zero)>=0,temp.abs(),BigInteger.valueOf(1));
                }
            gauss();
        }
    }
    public static void main(String args[])throws Exception {
        new Main();
    }
}
