/*
    扩展欧几里得算法
    计算给定区间内的不定方程ax+by+c=0的所有解的个数
    注意点：负数整除的进位问题，因此wa到爆了……
*/
#include <cstdio>
#include <algorithm>
#include <iostream>
using namespace std;
typedef long long LL;
LL ex_gcd(LL a, LL b, LL &x, LL &y) {
    if (b == 0) { x = 1, y = 0; return a; }
    LL r = ex_gcd(b, a%b, x, y);
    LL t = x;
    x = y; y = t - a/b *y;
    return r;
}
int main() {
    LL a, b, c, x1, x2, y1, y2;
    cin >> a>>b>>c>>x1>>x2>>y1>>y2;
    c = -c;
    if (a == 0 && b == 0) {
        if (c == 0) cout << (x2-x1+1)*(y2-y1+1) << endl;
        else cout << 0 << endl;
        return 0;
    }
    if (a == 0) {
        if (c % b == 0 && y1 <= c/b && c/b <= y2) cout << x2-x1+1 << endl;
        else cout << 0 << endl;
        return 0;
    }
    if (b == 0) {
        if (c % a == 0 && x1 <= c/a && c/a <= x2) cout << y2-y1+1 << endl;
        else cout << 0 << endl;
        return 0;
    }

    LL x, y, gcd;
    gcd = ex_gcd(a, b, x, y);
    if (c % gcd != 0) {
        cout << 0 << endl;
        return 0;
    }

    c /= gcd;
    x *= c; x1 -= x; x2 -= x;
    y *= c; y1 -= y; y2 -= y;
    c = -y2; y2 = -y1; y1 = c;
    b /= gcd;
    a /= gcd;

    if (b < 0) { c = -x2; x2 = -x1; x1 = c; b = -b; }
    if (a < 0) { c = -y2; y2 = -y1; y1 = c; a = -a; }
    LL t1, t2, r1, r2;
    t1 = (x1<0 || x1%b == 0) ? x1/b : x1/b+1;
    t2 = (x2>0 || x2%b == 0) ? x2/b : x2/b-1;
    r1 = (y1<0 || y1%a == 0) ? y1/a : y1/a+1;
    r2 = (y2>0 || y2%a == 0) ? y2/a : y2/a-1;

    cout << min(t2, r2)-max(t1, r1) + 1 << endl;

    return 0;
}
