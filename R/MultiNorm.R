####Norm on []X[]X[] space
MultiNorm<-function(v,p,tp,interval = c(0,1)){
	L = length(v)
	S = 0
	for(i in 1:p){
		ind = (1+(i-1)*(L/p)) :(i*(L/p))
		S = S + FakeInt(v[ind],tp,interval) 
	}
	S
}