library(reshape2)
library(ggplot2)
library(gridExtra)
reproduction = function(pop, p, q){
	genos=c()
	for(i in 1:pop){
		geno = sample(c("p", "q"), 2, rep=TRUE, prob=c(p,q))
		geno = paste(geno[1],geno[2], sep="");
		if(geno=="qp"){
			geno="pq"
		}
		genos = append(genos, geno)
	}
	genofreq = table(genos)
	pp = as.numeric(genofreq[1])/pop
	pq = as.numeric(genofreq[2])/pop
	qq = as.numeric(genofreq[3])/pop
	c(pp, pq, qq)
}		
evolution = function(pop, gen, p, q, wpp, wpq, wqq, pmut, qmut){
	i = 1;
	#Create first batch of offspring
	offspring = reproduction(pop, p, q)
	pp=offspring[1]; pq=offspring[2]; qq=offspring[3];
	allelehistory = cbind(i, p, q);
	genohistory = cbind(i, pp, pq, qq)
	while(i < gen){
		if((p+q)!=1){
			print(paste("Alleles do not add to 1!", p+q))
			q=1-p
		}
		#Select on Generation
		averagefitness = (pp*wpp)+(pq*wpq)+(qq*wqq)
		newpp = pp*wpp/averagefitness
		newpq = pq*wpq/averagefitness
		newqq = qq*wqq/averagefitness
		if(is.na(averagefitness)){
			print("er");
		}
		#Calculate allele frequency in gametes
		p=newpp+(.5*newpq); q=newqq+(.5*newpq);
		#Mutatate offspring
		#p = p-(p*pmut)+(q*qmut);
		#q = q+(p*pmut)-(q*qmut);
		#Create next generation
		offspring = reproduction(pop, p, q)
		pp=offspring[1]; pq=offspring[2]; qq=offspring[3];
		#Add to history
		allelehistory = rbind(allelehistory, cbind(i,p,q))
		genohistory = rbind(genohistory, cbind(i, pp, pq, qq))
		i=i+1;
	}
	colnames(genohistory) = c("generation", "pp", "pq", "qq")
	colnames(allelehistory) = c("generation", "p", "q")
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