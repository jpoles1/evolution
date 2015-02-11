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
#Used to calculate the equilibriam frequency at which a given deleterious allele is created and destroyed at equal rates by mutation and selection respectively.
#Selection coefficient is equal to 1-fitness
#equilibrium frequency is equal to the sqrt of mutation rate to lethal allele / selection coefficient
calculateMutationSelectionBalance = function(mutationRateLethal, fitness){
  sqrt(mutation/(1-fitness))
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
mutate = function(p, q, pmut, qmut){
	newp = p-(p*pmut)+(q*qmut);
	newq = q+(p*pmut)-(q*qmut);
	c(newp,newq)
}
evolution = function(gen, p, q, wpp, wpq, wqq, pmut, qmut){
	i = 1;
	allelehistory = cbind(i, p, q);
	genohistory = cbind(i, p^2, 2*p*q, q^2);
	while(i < gen){
		selectedalleles = selection(p, q, wpp, wpq, wqq);
		mutatedalleles = mutate(selectedalleles[1], selectedalleles[2], pmut, qmut)
		p=mutatedalleles[1]; q=mutatedalleles[2]
		allelehistory = rbind(allelehistory, cbind(i,p,q))
		genohistory = rbind(genohistory, cbind(i, p^2, 2*p*q, q^2))
		i=i+1;
	}
	colnames(allelehistory) = c("generation", "p", "q")
	colnames(genohistory) = c("generation", "pp", "pq", "qq")
	list(allelehistory, genohistory)
}
evolve = function(p = .5, gen = 50, wpp = 1, wpq = 1, wqq = 1, pmut=0, qmut=0){
	q=1-p
	x = evolution(gen, p, q, wpp, wpq, wqq, pmut, qmut)
	alleles = x[[1]]
	genotypes = x[[2]]
	allelemax = round(max(alleles[gen,2:3]), 2)
	allelemin = round(min(alleles[gen,2:3]), 2)
  alleledat = melt(alleles[,2:3])
  ggplot(alleledat, aes(x=alleledat$Var1, y=alleledat$value, col=alleledat$Var2))+ylim(0,1)+geom_point()+ggtitle(paste("Alleles Over Time; Final p=", round(alleles[gen,2],2), "q=", round(alleles[gen,3],2)))+xlab("Generation #")+ylab("Allele Freq.")+scale_color_discrete(name="Alleles", labels=c("p", "q"))
	plot2 = ggplot(melt(genotypes[,2:4]), aes(x=Var1, y=value, col=Var2))+ylim(0,1)+geom_point()+ggtitle(paste("Genotypes Over Time; Final pp=", round(genotypes[gen,2],2), "pq=", round(genotypes[gen,3],2), "qq=",round(genotypes[gen,4],2)))+xlab("Generation #")+ylab("Genotype Freq.")+scale_color_discrete(name="Genotypes", labels=c("pp", "pq", "qq"))
	grid.arrange(plot1, plot2, ncol=2, main=paste("Evolution with fitness: ","pp=",wpp,"; pq=",wpq,"; qq=",wqq, "\n and mutation rates p=", pmut, "; q=",qmut, sep=""))
}
changeFitness = function(){
	fit = 1
	while(fit>0){
		evolve(p=.5, gen = 50, wpp=fit)
		fit = fit-.05
		Sys.sleep(1.5)
	}
}