library(reshape2)
library(ggplot2)
library(gridExtra)
genotypeFromAllele = function(p, q){
	pp = p^2
	pq = 2*p*q
	qq = q^2
	c(pp, pq, qq)
}
alleleFromGenotype = function(pp, pq, qq){
	p = (pp) + (.5*pq);
	q = (qq) + (.5*pq);
	c(p, q)
}
selection = function(p, q, wpp, wpq, wqq){
	x = genotypeFromAllele(p,q)
	pp = x[1];pq = x[2];qq = x[3];
	averagefitness = (pp*wpp)+(pq*wpq)+(qq*wqq)
	newpp = pp*wpp/averagefitness
	newpq = pq*wpq/averagefitness
	newqq = qq*wqq/averagefitness
	y = alleleFromGenotype(newpp, newpq, newqq)
	p=y[1]; q=y[2];
	c(p,q)
}
evolution = function(gen, p, q, wpp, wpq, wqq){
	i = 1;
	allelehistory = cbind(i, p, q);
	genohistory = cbind(i, p^2, 2*p*q, q^2);
	while(i < gen){
		newalleles = selection(p, q, wpp, wpq, wqq)
		p=newalleles[1]; q=newalleles[2]
		allelehistory = rbind(allelehistory, cbind(i,p,q))
		genohistory = rbind(genohistory, cbind(i, p^2, 2*p*q, q^2))
		i=i+1;
	}
	colnames(allelehistory) = c("generation", "p", "q")
	colnames(genohistory) = c("generation", "pp", "pq", "qq")
	list(allelehistory, genohistory)
}
evolve = function(p = .5, gen = 50, wpp = 1, wpq = 1, wqq = 1){
	q=1-p
	x = evolution(gen, p, q, wpp, wpq, wqq)
	alleles = x[[1]]
	genotypes = x[[2]]
	#plot(alleles[,1], alleles[,2], col=2, ylim=c(0,1))
	#points(alleles[,3], col=3)
	plot1 = ggplot(melt(alleles[,2:3]), aes(x=Var1, y=value, col=Var2))+geom_point()
	plot2 = ggplot(melt(genotypes[,2:4]), aes(x=Var1, y=value, col=Var2))+geom_point()
	grid.arrange(plot1, plot2, ncol=2, main=paste("Evolution with fitness",wpp))
	#melt(alleles[,2:3])
}
fit = 1
while(fit>0){
	evolve(p=.5, gen = 50, wpp=fit)
	fit = fit-.05
	Sys.sleep(1.5)
}