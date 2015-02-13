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
	if((p+q)!=1){
		print(paste("error, doesn't add to 1", p, q))
	}
	c(p, q)
}
geneticdrift = function(pop, p, q){
	if((p+q)==1){
		
	}
	else{
		print(paste("Errors!", p, q))
	}
	c(p, q)
}
#Used to calculate the equilibriam frequency at which a given deleterious allele is created and destroyed at equal rates by mutation and selection respectively.
#Selection coefficient is equal to 1-fitness
#equilibrium frequency is equal to the sqrt of mutation rate to deleterious allele / selection coefficient
calculateMutationSelectionBalance = function(mutationRateLethal, fitness){
  sqrt(mutation/(1-fitness))
}
#Simulates migration using the fraction of the new population derived from migrant individuals
migrationByFraction = function(p, q, migrationrate){
  
}
reproduction = function(p, q, wpp, wpq, wqq, pmut, qmut, pop){
	#Get Genotypes of Parents
	x = genotypeFromAllele(p,q)
	#Select on Parents
	pp = x[1];pq = x[2];qq = x[3];
	averagefitness = (pp*wpp)+(pq*wpq)+(qq*wqq)
	newpp = pp*wpp/averagefitness
	newpq = pq*wpq/averagefitness
	newqq = qq*wqq/averagefitness
	#Calculate allele frequency in gametes
	y = alleleFromGenotype(newpp, newpq, newqq)
	p=y[1]; q=y[2];
	#Perform Mating
	alleledist = tryCatch({sample(c("p", "q"),pop,rep=TRUE,prob=c(p, q))}, 
		error = function(e){
			print(paste("Error", p, q))
		}
	)
	allelefreq = prop.table(table(alleledist))
	p = as.numeric(allelefreq[1])
	q = as.numeric(allelefreq[2])
	#Mutatate offspring
	p = p-(p*pmut)+(q*qmut);
	q = q+(p*pmut)-(q*qmut);
	c(p,q)
}
evolution = function(pop, gen, p, q, wpp, wpq, wqq, pmut, qmut){
	i = 1;
	allelehistory = cbind(i, p, q);
	genohistory = cbind(i, p^2, 2*p*q, q^2);
	while(i < gen){
		offspring = reproduction(p, q, wpp, wpq, wqq, pmut, qmut, pop);
		p=offspring[1]; q=offspring[2]
		allelehistory = rbind(allelehistory, cbind(i,p,q))
		genohistory = rbind(genohistory, cbind(i, p^2, 2*p*q, q^2))
		i=i+1;
	}
	colnames(allelehistory) = c("generation", "p", "q")
	colnames(genohistory) = c("generation", "pp", "pq", "qq")
	list(allelehistory, genohistory)
}
evolve = function(pop=1000, gen = 50, p = .5, wpp = 1, wpq = 1, wqq = 1, pmut=0, qmut=0){
	q=1-p
	x = evolution(pop, gen, p, q, wpp, wpq, wqq, pmut, qmut)
	alleles = x[[1]]
	genotypes = x[[2]]
	allelemax = round(max(alleles[gen,2:3]), 2)
	allelemin = round(min(alleles[gen,2:3]), 2)
  	#alleledat = melt(alleles[,2:3])
	plot1 = ggplot(melt(alleles[,2:3]), aes(x=Var1, y=value, col=Var2))+ylim(0,1)+geom_point()+ggtitle(paste("Alleles Over Time\n Final p=", round(alleles[gen,2],3), "; q=", round(alleles[gen,3],3), sep=""))+xlab("Generation #")+ylab("Allele Freq.")+scale_color_discrete(name="Alleles", labels=c("p", "q"))
	plot2 = ggplot(melt(genotypes[,2:4]), aes(x=Var1, y=value, col=Var2))+ylim(0,1)+geom_point()+ggtitle(paste("Genotypes Over Time\n Final pp=", round(genotypes[gen,2],3), "; pq=", round(genotypes[gen,3],3), "; qq=",round(genotypes[gen,4],3), sep=""))+xlab("Generation #")+ylab("Genotype Freq.")+scale_color_discrete(name="Genotypes", labels=c("pp", "pq", "qq"))
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