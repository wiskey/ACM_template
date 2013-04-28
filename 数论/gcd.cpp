// gcd:最大公约数
int gcd(int a, int b) {
	if (b == 0) return a;
	return gcd(b, a%b);
}

// extern_gcd 扩展欧几里得算法。解ax + by = gcd(a, b)的一组解x，y
//返回值是a，b的最大公约数，如果不需要求gcd则可以不返回。
int extern_gcd(int a, int b, int &x, int &y) {
	if (b == 0) { x = 1, y = 0; return a; }
	int r = extern_gcd(b, a%b, x, y);
	int t = x;
	x = y; y = t - a/b*y;
	return r;
}