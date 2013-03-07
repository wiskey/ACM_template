#include <cstdio>

//Euler函数：对于一个正整数n，小于n且和n互质的正整数（包括 1）的个数，记作 φ(n)。
//欧拉定理：对于互质的正整数 a 和 n ，有 a^φ(n) ≡ 1 mod n 

int Euler(int n) {
    int e = n;
    for (int i=2; i*i<=n; i++)
        if (n % i == 0) {
            e = e*(i-1)/i;
            while (n % i == 0) n /= i;
        }
    if (n > 1) e = e*(n-1) / n;
    return e;
}
int main() {
    int n;
    while (scanf("%d", &n) == 1)
        printf("%d\n", Euler(n));
    return 0;
}
