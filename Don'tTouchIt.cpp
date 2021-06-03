#include<bits/stdc++.h>
using namespace std ;

#define ll long long
#define pb push_back
#define all(v) v.begin(),v.end()
#define ff first
#define ss second
#define ar array
#define endl '\n'
const int mod = 1e9+7 ;

int top_46prime[46] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
					    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
						109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
						173, 179, 181, 191, 193, 197, 199 } ;

/*        | 1                                  if n == 0
	a^n = | ( a ^ ( n/2 ) ) ^ 2                if n is even & n > 0
		  | ( ( a ^ ( (n-1) / 2 ) ^ 2 ) * a    if n is odd  & n > 0
*/
ll binpow(ll a , ll b ){
	a %= mod ;
	ll res = 1 ;
	while( b > 0 ) {
		if( b & 1 )
			res = res * a % mod ;
		a = a*a % mod ;
		b >>= 1;
	}
	return res ;
}
/* gcd( a , b ) = | a           if b == 0
				  | gcd(b,a%b)  if b > 0
*/

int gcd ( int a , int b ){
	while( b ){
		a %= b ;
		swap(a,b) ;
	}
	return a ;
}

//  ( a * x ) + ( b * x ) = gcd( a , b )

int extgcd( int a , int b , int& x , int& y ){
	if( b == 0 ) {
		 x = 1 ;
		 y = 0 ;
		 return a ;
	}

	int x1 , y1 ;
	int d = extgcd(b , a%b , x1 , y1) ;
	x = y1 ;
	y = x1 - y1 * ( a / b ) ;
	return d ;
}

// ----------------------------Sieve of Eratosthenes here
const int MxN = (int)1e6 ;
int lp[MxN+1] ;
vector<int> primes ;

void cal_primes(){
	for( int i = 2 ; i <= MxN ; ++i ){
		 if( lp[i] == 0 ) {
			 lp[i] = i ;
			 primes.push_back(i) ;
		 }
		 for( int j = 0 ; j < (int)primes.size() && primes[j] <= lp[i] && i*primes[j] <= MxN ; ++j ){
			  lp[i * primes[j]] = primes[j] ;
		 }
	}
}
// ------------------------------factorial here

ll f1[MxN+1], f2[MxN+2], inv[MxN+1] ;

void cal_fact(){
	f1[0] = f2[0] = 1 ;
	inv[1] = 1 ;
	for( int i = 2 ; i <= MxN ; ++i)
		inv[i] = mod - (mod/i) * inv[mod%i] % mod ;

	for( int i = 1 ; i <= MxN; i++){
		f1[i] = f1[i-1] * i % mod ;
		f2[i] = f2[i-1] * inv[i] % mod ;
	}
	// NcK => f1[N] * f2[K] % mod * f2[N-k] % mod
	
}

pair<ll, ll> fib(ll n){
	if ( n == 0 )
	   return {0, 1} ;

	auto p = fib( n >> 1 ) ;
	ll c = p.first * ( 2 * p.second - p.first ) % mod ;
	ll d = ( p.first * p.first % mod + p.second * p.second % mod ) % mod ;
	if( n & 1 )
		return {d, c+d} ;
	else
		return {c, d} ;
}
// ----------- Bit Manipulation --------------

// set the i'th bit to 1
int setBit( int set , int i ){
	return set | ( 1 << i ) ;
}

// set the i'th bit to 0
int unsetBit( int s , int i ){
	return s & ~( 1 << i ) ;
}

// Checks if the i'th bit is set
bool isSet( int s , int i ){
	return ( s & ( 1 << i ) ) != 0 ;
}

// Toggle the i'th bit to zero
int toggleBit( int s , int i ){
	return s ^ ( 1 << i ) ;
}

// return a number with the first n bit set to 1
int setAll(int n){
	return ( 1 << n ) - 1 ;
}

// verifies if a number n is power of 2
bool isPowerOfTwo(int n){
	return n > 0 && ( n & (n-1) ) == 0 ;
}

/*--------- Z-Algorithm ---------*/
vector<int> calculateZ(string x, string y){
	
	string z = x + "#" + y ;
	int n = z.length() ;
	vector<int> Z(n) ;
	int l = 0 , r = 0 , k ;

	for( int i = 1 ; i < n ; i++){
		if( i > r ){
			l = r = i ;
			while( r < n && z[r-l] == z[r] ){
				r++ ;
			}
			Z[i] = r - l ;
			r-- ;
		}
		else{
			k = i - l ;
			if( Z[k] < r - i + 1 ){
				Z[i] = Z[k] ;
			}
			else{
				l = i ;
				while( r < n && z[r-l] == z[r]){
					r++ ;
				}
				Z[i] = r-l ;
				r-- ;
			}
		}
	}
	return Z ;
}

/*------------------------Segment Tree------------------------

struct segtree {

	int size ;
	vector<ll> seg ;

	void init(int n){
		size = 1 ;
		while( size < n ){
			size *= 2 ;
		}
		seg.assign(2*size,0ll) ;
	}
	void build(vector<int> &a){
		build(a,0,0,size) ;
	}

	void build(vector<int> &a, int x, int lx, int rx){
		if( rx - lx == 1 ){
			if( lx < (int)a.size() ){
				seg[x] = a[lx] ;
			}
			return ;
		}
		{
			int m = (rx+lx)/2 ;
			build(a, 2*x+1, lx, m ) ;
			build(a, 2*x+2, m, rx) ;
			seg[x] = seg[2*x+1] + seg[2*x+2] ;
		}
	}

	void update(int i, int val){
		update(i, val, 0, 0, size) ;
	}

	void update(int i, int val, int x, int lx, int rx){
		if( rx - lx == 1 ){
			seg[x] = val ;
		}
		else {
			int m = (rx+lx)/2 ;
			if( i < m )
				update(i, val ,2*x+1, lx, m ) ;
			else
				update(i, val, 2*x+2, m, rx) ;
			seg[x] = seg[2*x+1] + seg[2*x+2] ;
		}
	}

	ll query(int a, int b){
		return query(a,b,0,0,size) ;
	}

	ll query(int a, int b, int x, int lx, int rx){
		if( a >= rx || b <= lx )
			return 0 ;
		if( a <= lx && b >= rx )
			return seg[x] ;

		int m = (lx+rx)/2 ;
		ll s1 = query(a,b,2*x+1,lx,m) ;
		ll s2 = query(a,b,2*x+2,m,rx) ;
		return s1+s2 ;
	}

};
//------------------------part2

const int Mxn = 2e5 ;
int n, q ;
ll a[Mxn] ;
struct segtree{
	ll x ;
} st[1<<19] ;
void upd( int l1, ll val, int i = 1, int l2 = 0, int r2 = n-1){
	if( l2 == r2 ){
		st[i].x = val ;
		return ;
	}
	int m2 = (l2+r2)/2 ;
	if( l1 <= m2 ){
		upd( l1, val, 2*i, l2, m2 ) ;
	}
	else{
		upd( l1, val, 2*i+1, m2+1, r2) ;
	}
	st[i].x = (st[2*i].x^st[2*i+1].x) ;
}
ll qry( int l1, int r1, int i = 1, int l2 = 0 , int r2 = n-1 ){
	if( l1 <= l2 && r1 >= r2 ){
		return st[i].x ;
	}
	int m2 = ( r2 + l2 ) / 2 ;
	ll l = 0 ;
	ll r = 0;
	if( l1 <= m2 ){
		l = qry(l1, r1, 2*i, l2, m2 ) ;
	}
	if( r1 > m2 )
		r = qry( l1, r1, 2*i+1, m2+1, r2) ;
	return  l^r ;
}
------------------------------------------------------------------------------------*/

// sum of squares = (2*n+1)(n+1)(n) / 6


string laxi_max_str_after_rotation_right_to_left(string s){
		s += s ;
        int i = 0, j = 1, k = 0 ;
        int n = s.size() ;
        while (j + k < n){
            if (s[i+k] == s[j+k]){
                k += 1 ;
                continue ;
			}
            else if (s[i+k] > s[j+k]){
                j = j + k + 1 ;
			}
            else{
                i = max(i + k + 1, j) ;
                j = i + 1 ;
			}
            k = 0;
        }
        return s.substr(i)  ;
}


int dx[8] = { 0, 1, 0, -1, 1, 1, -1, -1 };
int dy[8] = { 1, 0, -1, 0, 1, -1, 1, -1 };

void solve(){

}


int main(){
	int t ; cin >> t ;

	//clock_t start , finish ;
	//start = clock() ;
	while(t--){
		solve() ;
	}
	//finish = clock() ;
	//cout << fixed << setprecision(3) << (float)(finish - start)/CLOCKS_PER_SEC ;
	
}

