#include <stdio.h>
#include <string.h>

#define MAXN 1005
#define MOD 1000000007
#define INV2 500000004

long long dp[MAXN][MAXN];

int main()
{
    int T, t, n, s, i, j;
    scanf("%d", &T);
    for(t=1; t<=T; t++)
    {
        scanf("%d%d", &n, &s);
        memset(dp, 0, sizeof dp);
        for(i=2, dp[1][1]=1; i<=n && scanf("%d", &s); i++)
            if(s && i!=n)
                for(j=1; j<=i; j++)
                    dp[i][j] = dp[i-1][j-1];
            else
            {
                for(j=i; j; j--)
                    dp[i][j] = (dp[i][j+1] + dp[i-1][j]) * INV2 % MOD;
                dp[i][1] <<= 1;
            }
        printf("Case #%d: %I64d\n", t, dp[n][1] * INV2 % MOD);
    }
    return 0;
}
