int lengthOfLIS(vector<int>& a) {
        vector<int> dp ; 
        for( auto x : a ){
            auto it = lower_bound(dp.begin(),dp.end(),x) ; 
            if(it == dp.end() ){
                dp.push_back(x) ;
            }   
            else{
                int in = it-dp.begin() ; 
                dp[in] = x ;
            }
        }
        return dp.size() ;
    }