/*
    裸最大公约数
*/
#include <cstdio>
#include <cstring>
using namespace std;
int gcd(int a, int b) {
    if (b == 0) return a;
    return gcd(b, a%b);
}
int main() {
    int n, c = 1;
    scanf("%d", &n);
    for (int i=2; i<n; i++) if (gcd(i, n) == 1) c++;
    printf("%d\n", c);

    return 0;
}
