#include <stdio.h>
#include <algorithm>

using namespace std;

#define NTT_MOD 998244353	//ģ��998244353=7*17*2^23+1
#define NTT_G 3			//��ģ����Ӧ��ԭ��ֵ
#define NTT_MAXW 30		//ԭ���������󳤶�,ȡlog2(MOD)����

inline int ntt_pow(long long a, int n = NTT_MOD - 2)
{	//������,Ĭ������Ԫ
	long long res = 1;
	do
	{
		if (n & 1)
			res = res * a % NTT_MOD;
		a = a * a % NTT_MOD;
	} while (n >>= 1);
	return (int)res;
}

void ntt(int* x, int n, bool inv = false)
{	//�������۱任,inv��ʾ�Ƿ���з��任,Ĭ��false��ʾ���任,�任����ӦΪ2�ı���,�����Ծ����������г�����
	int i, j, k, d, w, cur, tmp;
	static int g[NTT_MAXW], ng[NTT_MAXW];	//ԭ��������
	if (!*ng)	//��ʼ��
		for (i = 1, k = NTT_MOD - 1 >> 1, *ng = ntt_pow(NTT_G); !(k & 1); i++, k >>= 1)
		{
			g[i] = ntt_pow(NTT_G, k);
			ng[i] = ntt_pow(*ng, k);
		}
	for (i = n >> 1, j = 1; j < n; j++, i ^= k)
	{	//����λ��
		if (i < j)
			swap(x[i], x[j]);
		for (k = n >> 1; i&k; k >>= 1)
			i ^= k;
	}
	for (d = 2, k = 1; d <= n; d <<= 1, k++)	//��������
		for (i = 0,w = inv ? ng[k] : g[k]; i < n; i += d)
			for (j = i, cur = 1; j < i + (d >> 1); j++)
			{
				tmp = (long long)x[j + (d >> 1)] * cur % NTT_MOD;
				x[j + (d >> 1)] = (x[j] - tmp + NTT_MOD) % NTT_MOD;
				x[j] = (x[j] + tmp) % NTT_MOD;
				cur = (long long)cur * w % NTT_MOD;
			}
	if (inv)	//����任�ĳ�����Ԫ����
		for (i = 0, w = ntt_pow(n); i < n; i++)
			x[i] = (long long)x[i] * w % NTT_MOD;
}

#define MAXN 100005

int a[MAXN << 2];
int b[MAXN << 2];
long long ans[MAXN];
long long jc[MAXN] = { 1 };
long long ic[MAXN] = { 1 };
long long p2[MAXN] = { 1 };

int main()
{
	int t, n, i, len;
	scanf("%d", &t);
	for (i = 1; i<MAXN; i++)
	{
		jc[i] = jc[i - 1] * i % NTT_MOD;
		ic[i] = ntt_pow(jc[i], NTT_MOD - 2);
		p2[i] = p2[i - 1] * 2 % NTT_MOD;
	}
	while (t-- && scanf("%d", &n) > 0)
	{
		for (i = 0; i<n; i++)
			scanf("%d", a + i);
		sort(a, a + n);
		for (len = 1; len<n << 1; len <<= 1);
		for (i = 0; i<n; i++)
		{
			a[i] = p2[i] * a[i] % NTT_MOD * jc[n - 1 - i] % NTT_MOD;
			b[i] = ic[i];
		}
		do a[i] = b[i] = 0;
		while (++i < len);
		ntt(a, len, 0);
		ntt(b, len, 0);
		for (i = 0; i<len; i++)
			b[i] = (long long)a[i] * b[i] % NTT_MOD;
		ntt(b, len, true);
		for (i = n, ans[n] = 0; i--;)
			printf("%d ", ans[i] = (ans[i + 1] + (long long)b[i] * ic[n - 1 - i] % NTT_MOD) % NTT_MOD);
		putchar('\n');
	}
	return 0;
}