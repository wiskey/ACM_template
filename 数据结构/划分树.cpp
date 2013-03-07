/*
    带sum域维护的划分树
    HDU 3473 Minimum Sum
*/
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;

const int MAXN = 100010;
int sorted[MAXN]={0};   //对原集合中元素排序后的值
int val[20][MAXN]={0};  //val记录第k层当前位置的值
int toleft[20][MAXN]={0}; //记录元素所在区间当前位置前的元素进入到左子树的个数
int lnum, rnum;         //询问区间里面k-th数左侧和右侧的数的个数
long long sum[20][MAXN]={0};    //记录比当前元素小的元素的和
long long lsum, rsum;   //询问区间里面k-th数左侧数之和与右侧数之和
int n;

void build(int l, int r, int d) {
    if (l == r) return ;
    int mid = (l + r) >> 1;
    int same = mid - l + 1;
    for (int i=l; i<=r; i++)
        if (val[d][i] < sorted[mid])
            same--;
    int lp = l, rp = mid+1;
    for (int i=l; i<=r; i++) {
        if (i == l) {
            toleft[d][i] = 0;
            sum[d][i] = 0;
        } else {
            toleft[d][i] = toleft[d][i-1];
            sum[d][i] = sum[d][i-1];
        }

        if (val[d][i] < sorted[mid]) {
            toleft[d][i]++;
            sum[d][i] += val[d][i];
            val[d+1][lp++] = val[d][i];
        } else if (val[d][i] > sorted[mid])
            val[d+1][rp++] = val[d][i];
        else {
            if (same) {
                same--;
                toleft[d][i]++;
                sum[d][i] += val[d][i];
                val[d+1][lp++] = val[d][i];
            } else val[d+1][rp++] = val[d][i];
        }
    }
    build(l, mid, d+1);
    build(mid+1, r, d+1);
}
int query(int a, int b, int k, int l, int r, int d) {
    if (a == b) return val[d][a];

    int mid = (l + r) >> 1;
    int s, ss;
    long long sss;
    if (a == l) {
        s = toleft[d][b];
        ss = 0;
        sss = sum[d][b];
    } else {
        s = toleft[d][b] - toleft[d][a-1];
        ss = toleft[d][a-1];
        sss = sum[d][b] - sum[d][a-1];
    }

    if (s >= k) {
        a = l + ss;
        b = l + ss + s - 1;
        return query(a, b, k, l, mid, d+1);
    } else {
        lnum += s;
        lsum += sss;
        a = mid+1 + a - l - ss;
        b = mid+1 + b - l - toleft[d][b];
        return query(a, b, k-s, mid+1, r, d+1);
    }
}
int main() {
    freopen("in.txt", "r", stdin);
    int cas, n, m;
    long long s[MAXN]={0};

    scanf("%d", &cas);
    for (int t=1; t<=cas; t++) {
        scanf("%d", &n);
        s[0] = 0;
        for (int i=1; i<=n; i++) {
            scanf("%d", &sorted[i]);
            val[0][i] = sorted[i];
            s[i] = s[i-1] + sorted[i];
        }
        sort(sorted+1, sorted+1+n);
        build(1, n, 0);

        printf("Case #%d:\n", t);
        int a, b, k;
        scanf("%d", &m);
        while (m--) {
            scanf("%d %d", &a, &b);
            a++, b++;
            k = (b-a+2)>>1;
            lsum = lnum = 0;
            int ave = query(a, b, k, 1, n, 0);
            rnum = b-a+1 - lnum;
            rsum = s[b] - s[a-1] - lsum;
            printf("%I64d\n", rsum-lsum+(long long)(lnum-rnum)*(long long)ave);
        }
        printf("\n");
    }
    return 0;
}
