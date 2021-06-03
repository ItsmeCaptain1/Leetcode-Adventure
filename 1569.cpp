#define ll long long 
const ll mod = 1e9+7 ; 
const ll Mxn = 1e3+1 ;
bool flag = true ;

ll f1[Mxn+1], f2[Mxn+2], inv[Mxn+1] ; 

void cal(){
    f1[0] = f2[0] = 1 ; 
    inv[1] = 1 ;
	for( int i = 2 ; i <= Mxn ; ++i)
		inv[i] = mod - (mod/i) * inv[mod%i] % mod ;

	for( int i = 1 ; i <= Mxn; i++){
		f1[i] = f1[i-1] * i % mod ;
		f2[i] = f2[i-1] * inv[i] % mod ;
	}
}
class Solution {
public:
    ll go(int in, vector<int>& a, int l, int r){
        // cout << in << " " << l << " " << r << endl;
        if( l>=r ) return 1 ; 
        ll n = 0, m = 0, fl = 0, fh = 0 ;
        for( int i = in+1 ; i < a.size() ; i++ ){
            if( a[i] >= l && a[i] < a[in] ){
                n++ ;
                if( fl == 0 ) fl = i ;
            } 
            if( a[i] <= r && a[i] > a[in] ){ 
                m++ ;
                if( fh == 0 ) fh = i ;
            }
        }
        
        ll curr = (f1[n+m] * f2[n]%mod * f2[m]%mod)  ;
        cout << n << " " << m << endl ;
        cout << curr << endl ;        
        return (curr * go(fl,a,l,a[in]-1)) % mod * go(fh,a,a[in]+1,r) % mod ;
    }
    
    int numOfWays(vector<int>& nums) {
        
        if(flag){
            flag = true ;
            cal() ;
        }
        return go(0,nums,1,nums.size())-1 ;
    }
};