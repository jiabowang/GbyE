`GAPIT.Manhattan` <-
function(GI.MP = NULL,GD=NULL,name.of.trait = "Trait",plot.type = "Genomewise",Nenviron=2,Name_environ,
DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,interQTN=NULL,plot.style="Oceanic",CG=NULL,plot.bin=10^9){
    #Object: Make a Manhattan Plot
    #Options for plot.type = "Separate_Graph_for_Each_Chromosome" and "Same_Graph_for_Each_Chromosome"
    #Output: A pdf of the Manhattan Plot
    #Authors: Alex Lipka, Zhiwu Zhang, Meng Li and Jiabo Wang
    # Last update: Oct 10, 2016
	#Add r2 between candidata SNP and other markers in on choromosome
    ##############################################################################################
    #print("Manhattan ploting...")
    #Name_environ=c("E1","E2","E3","E4")
    #GI.MP=result[,-1]
    #name.of.trait="test"
    #seqQTN=NULL
    #interQTN=NULL
    #Nenviron=4
    ##cutOff=0.01
    #band=5
    #plot.style="Oceanic"
    #DPP=50000
    #print(seqQTN)
    #do nothing if null input
    if(is.null(GI.MP)) return
	if(length(Name_environ)!=Nenviron) print("@@@")
    position.only=T
    if(!is.null(seqQTN)|!is.null(interQTN)){
      
        #seqQTN=-seqQTN
        position.only=F
      
      
    }
    
    #if(is.null(GD)) print ("GD is not same dim as GM")
    total_chromo=max(GI.MP[,1])
    if(total_chromo%%(Nenviron)!=0) print("!!!")
    num_add_chromo=total_chromo/(Nenviron)
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
	GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    if(!is.null(interQTN))GI.MP[interQTN,borrowSlot]=2
    
    GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))
    GI.MP=GI.MP[order(GI.MP[,2]),]
    GI.MP=GI.MP[order(GI.MP[,1]),]
    
    
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    
    #Replace P the -log10 of the P-values
    
    GI.MP[,3] <-  -log10(GI.MP[,3])
    
    
    GI.MP[,5]=1:(nrow(GI.MP))
    y.lim <- ceiling(max(GI.MP[,3]))
    chm.to.analyze <- unique(GI.MP[,1])
    
    
    chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    add_chromo=chm.to.analyze[1:num_add_chromo]
    numCHR= num_add_chromo
    
    #Genomewise plot
    #if(plot.type == "Genomewise")
    #{
        #print("Manhattan ploting Genomewise")
        #Set corlos for chromosomes
        #nchr=max(chm.to.analyze)
        nchr=numCHR

    #Set color schem            
        ncycle=ceiling(nchr/band)
        ncolor=band*ncycle
        #palette(rainbow(ncolor+1))
        cycle1=seq(1,nchr,by= ncycle)
        thecolor=cycle1
        for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
      	col.Rainbow=rainbow(ncolor+1)[thecolor]     	
     	  col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
    	  col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))   	
		    col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
 		    col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2)) 		
 		    col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2)) 		
 		    col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
 		    #col.Oceanic=rep(c(	'#EC5f67',	'#F99157',	'#FAC863',	'#99C794',	'#5FB3B3',	'#6699CC',	'#C594C5',	'#AB7967'),ceiling(numCHR/8))
 		    #col.Oceanic=rep(c(	'#EC5f67',		'#FAC863',	'#99C794',		'#6699CC',	'#C594C5',	'#AB7967'),ceiling(numCHR/6))
 		    col.Oceanic=rep(c(	'#EC5f67',		'#FAC863',	'#99C794',		'#6699CC',	'#C594C5'),ceiling(numCHR/5))
 		    col.cougars=rep(c(	'#990000',		'dimgray'),ceiling(numCHR/2))
 		
        if(plot.style=="Rainbow")plot.color= col.Rainbow
        if(plot.style =="FarmCPU")plot.color= col.FarmCPU
        if(plot.style =="Rushville")plot.color= col.Rushville
        if(plot.style =="Congress")plot.color= col.Congress
        if(plot.style =="Ocean")plot.color= col.Ocean
        if(plot.style =="PLINK")plot.color= col.PLINK
 		    if(plot.style =="Beach")plot.color= col.Beach
 		    if(plot.style =="Oceanic")plot.color= col.Oceanic
 		    if(plot.style =="cougars")plot.color= col.cougars
 		
		#FarmCPU uses filled dots
    	mypch=1
    	if(plot.style =="FarmCPU")mypch=20
    	        
        GI.MP <- GI.MP[order(GI.MP[,2]),]
        GI.MP <- GI.MP[order(GI.MP[,1]),]

        GI.MP[,6]=1:(nrow(GI.MP))
        for(i in 1:(Nenviron))
        {
          GI.MP[(1+(i-1)*(nrow(GI.MP)/Nenviron)):((nrow(GI.MP)/Nenviron+(i-1)*(nrow(GI.MP)/Nenviron))),6]=i
        }
        
            #result[(1+(i-1)*(nrow(result)/n_e)):((nrow(result)/n_e+(i-1)*(nrow(result)/n_e))),2]=chrom

        
        
        #print("Manhattan data sorted")
        #print(chm.to.analyze)
        
        #change base position to accumulatives (ticks)
        #plot by number of environments
        add_store=GI.MP[GI.MP[,6]==1,]

        x_store=NULL
        z_store=NULL
       for(k in 1:(Nenviron))
    {
     
        MP_store=GI.MP[GI.MP[,6]==k,]
        
        MP_store[,1]=add_store[,1]
        MP_store[,2]=add_store[,2]
        index_GI=MP_store[,3]>0
        MP_store <- MP_store[index_GI,]
        ticks=NULL
        lastbase=0
        for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
        #print("Manhattan chr processed")
        #print(length(index))
        #print(length(ticks))
        #print((ticks))
        #print((lastbase))
        
        x0 <- as.numeric(MP_store[,2])
        y0 <- as.numeric(MP_store[,3])
        z0 <- as.numeric(MP_store[,1])
        position=order(y0,decreasing = TRUE)
        #index0=GAPIT.Pruning(y0[position],DPP=DPP)
        values=y0[position]
        if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
            }else{
         
        
        values=sqrt(values)  #This shift the weight a little bit to the low building.

        #Handler of bias plot
        rv=runif(length(values))
        values=values+rv
        values=values[order(values,decreasing = T)]

        theMin=min(values)
        theMax=max(values)
        range=theMax-theMin
        interval=range/DPP

        ladder=round(values/interval)
        ladder2=c(ladder[-1],0)
        keep=ladder-ladder2
        #index=which(keep>0)


        index=position[which(keep>0)]
            }
        
        
        x=x0[index]
        y=y0[index]
        z=z0[index]
        #Draw circles with same size and different thikness
        
        
        #Label 
        index_bon=y>=bonferroniCutOff
        x_store=append(x_store,x[index_bon])
        z_store=append(z_store,z[index_bon])
        
        #print("Manhattan done Genomewise")
    }#Nenvironments
        pdf(paste("GAPIT.", name.of.trait,".Manhattan.Plot.Genomewise.pdf" ,sep = ""), width = 13,height=15)
        par(mfrow=c(Nenviron,1))
        #par(mar = c(0.1,6,0.1,6))
        
    for(k in 1:(Nenviron))
    {
        if(k==Nenviron){#par(mfrow=c(Nenviron,1))
        par(mar = c(3,6,0.5,6))
        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(0,6,0.5,6))
        
        }
     
        MP_store=GI.MP[GI.MP[,6]==k,]
        #print(k)
        #print(head(MP_store))
        #print(head(add_store))
        MP_store[,1]=add_store[,1]
        MP_store[,2]=add_store[,2]
        MP_store[,4]=add_store[,4]
        
        index_GI=MP_store[,3]>0
        MP_store <- MP_store[index_GI,]
        ticks=NULL
        lastbase=0
        for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
        #print("Manhattan chr processed")
        #print(length(index))
        #print(length(ticks))
        #print((ticks))
        #print((lastbase))
        
        x0 <- as.numeric(MP_store[,2])
        y0 <- as.numeric(MP_store[,3])
        z0 <- as.numeric(MP_store[,1])
        position=order(y0,decreasing = TRUE)
        #index0=GAPIT.Pruning(y0[position],DPP=DPP)
        values=y0[position]
        if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
            }else{
         
        
        values=sqrt(values)  #This shift the weight a little bit to the low building.

        #Handler of bias plot
        rv=runif(length(values))
        values=values+rv
        values=values[order(values,decreasing = T)]

        theMin=min(values)
        theMax=max(values)
        range=theMax-theMin
        interval=range/DPP

        ladder=round(values/interval)
        ladder2=c(ladder[-1],0)
        keep=ladder-ladder2
        #index=which(keep>0)


        index=position[which(keep>0)]
            }
        
        
        x=x0[index]
        y=y0[index]
        z=z0[index]

        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
        QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        interQTN=MP_store[which(MP_store[,borrowSlot]==2),]
        #Draw circles with same size and different thikness
        size=1 #1
        ratio=10 #5
        base=1 #1
        if(k==1)themax=ceiling(max(y))
        if(k!=1)
        {
            if(themax<ceiling(max(y)))themax=ceiling(max(y))
        }
        themin=floor(min(y))
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        s=size-wd/ratio/2
        
        #print("Manhattan XY created")
       ####xiaolei update on 2016/01/09 
            
        	plot(y~x,xlab="",ylab=expression(-log[10](italic(p))) ,ylim=c(0,themax),
        	cex.axis=1.5, cex.lab=2, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=wd,cex=s+.3,cex.main=2.5)
            
        #Label QTN positions
        #print(head(QTN))
        #print(head(interQTN))
          if(position.only){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
            #print("$$$$$$")
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(interQTN[,2], interQTN[,3], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        
        #}
        abline(v=x_store,col=plot.color[z_store],lty=3,untf=T,lwd=1.9)

        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,col="forestgreen")
        #mtext("n= 800",side=4,line=2.5,col="black")
        #Set axises
        #print(themax)
        axis(2, at=1:themax,cex.axis=1.5,labels=1:themax,tick=F)
        if(k==Nenviron) axis(1, at=ticks,cex.axis=1.5,labels=chm.to.analyze,tick=F)
        
        if(k==1)mtext(side=4,paste("ADD",name.of.trait,sep="  "),line=3)
        if(k!=1)mtext(side=4,paste("G&E ",Name_environ[k],sep="  "),line=3)
        box()
        palette("default")
        
        #print("Manhattan done Genomewise")
    }#Nenvironments
        #print(x[index_bon])
    #} #Genomewise plot
    dev.off()
    print("GAPIT.Manhattan accomplished successfully!zw")
} #end of GAPIT.Manhattan
#=============================================================================================
