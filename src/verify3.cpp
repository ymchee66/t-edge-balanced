// Independent check that an edge list is a 3-edge-balanced graph in G(n,k).
// Counts c (K3), p (P4), s (K_{1,3}), d (P3 u K2) from their definitions,
// m (3K2) = C(k,3) - c - p - s - d, and compares with the target profile
// c* = 4k(k-1)(k-2)/(3(n+1)(n^2-n-4)), p* = 3(n-3)c*, s* = (n-3)c*,
// d* = 3(n-3)(n-4)c*/2, m* = (n-3)(n-4)(n-5)c*/8.
// Usage: verify3 n k graph6file
#include <bits/stdc++.h>
using namespace std;
typedef __int128 L;
string str(L x){ if(x==0) return "0"; bool neg=x<0; if(neg) x=-x; string s; while(x){s+=char('0'+x%10);x/=10;} if(neg)s+='-'; reverse(s.begin(),s.end()); return s; }
int main(int argc,char**argv){
    int n=atoi(argv[1]); long long k=atoll(argv[2]);
    ifstream f(argv[3]); string g; f >> g;
    if(g.empty()){puts("BAD file");return 1;}
    size_t pos=0; int gn;
    if((unsigned char)g[0]==126){ if(g.size()<4){puts("BAD header");return 1;} gn=((g[1]-63)<<12)|((g[2]-63)<<6)|(g[3]-63); pos=4; }
    else { gn=g[0]-63; pos=1; }
    if(gn!=n){puts("BAD vertex count");return 1;}
    vector<pair<int,int>> E; long long bit=0;
    for(int j=1;j<n;j++) for(int i=0;i<j;i++,bit++){
        size_t byte=pos+bit/6; if(byte>=g.size()){puts("BAD length");return 1;}
        int val=g[byte]-63; if((val>>(5-bit%6))&1) E.push_back({i,j});
    }
    if((long long)E.size()!=k){puts("BAD edge count");return 1;}
    vector<vector<char>> A(n,vector<char>(n,0)); vector<vector<int>> nb(n); vector<long long> d(n,0);
    for(auto&e:E){ if(A[e.first][e.second]){puts("BAD repeated edge");return 1;} A[e.first][e.second]=A[e.second][e.first]=1; nb[e.first].push_back(e.second); nb[e.second].push_back(e.first); d[e.first]++; d[e.second]++; }
    L c=0; for(auto&e:E){ for(int w:nb[e.first]) if(A[e.second][w]) c++; } c/=3;
    L s=0; for(int x=0;x<n;x++) s+=(L)d[x]*(d[x]-1)*(d[x]-2)/6;
    L p=0; for(auto&e:E) p+=(L)(d[e.first]-1)*(d[e.second]-1); p-=3*c;
    L dd=0;  // P3 u K2: for each path x-y-z, count edges disjoint from {x,y,z}
    for(int y=0;y<n;y++){ auto&N=nb[y]; for(size_t i=0;i<N.size();i++) for(size_t j=i+1;j<N.size();j++){ int x=N[i],z=N[j];
        L touching=d[x]+d[y]+d[z]-2-(A[x][z]?1:0); dd+=(L)k-touching; } }
    L K3=(L)k*(k-1)*(k-2)/6; L m=K3-c-p-s-dd;
    L num=(L)4*k*(k-1)*(k-2), den=(L)3*(n+1)*((L)n*n-n-4);
    if(num%den){puts("target not integral");return 1;}
    L cs=num/den, ps=3*(L)(n-3)*cs, ss=(L)(n-3)*cs, ds=3*(L)(n-3)*(n-4)*cs/2, ms=(L)(n-3)*(n-4)*(n-5)*cs/8;
    bool ok=(c==cs&&p==ps&&s==ss&&dd==ds&&m==ms);
    printf("%d %lld %s c=%s p=%s s=%s d=%s m=%s\n",n,k,ok?"3-edge-balanced":"FAIL",str(c).c_str(),str(p).c_str(),str(s).c_str(),str(dd).c_str(),str(m).c_str());
    return ok?0:1;
}
