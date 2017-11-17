`G_E.Simulation` <-
function(GD,GM=NULL,ha2,hi2,he2,NE,NQTN,NETN,effectunit=1,indepent=TRUE){
    

    #############################
    #GD
    
    #############################
    if(ha2+hi2+he2>=1) 
    {
        print("simulation paremeter error in h2")
        exit
    }

    Y=matrix(0,nrow(GD),NE)
    X=GD[,-1]
    taxa=as.character(GD[,1])
    m=ncol(X)
    n=nrow(X)
    #duplicate taxa

    ####### simulate Addtive
    addeffect=rnorm(NQTN,0,1)
    QTN.position=sample(1:m,NQTN,replace=F)
    SNPQ=as.matrix(X[,(QTN.position)])
    Qeffect=SNPQ%*%addeffect
    Qeffectvar=var(Qeffect)
    
    ####### calculate P A I E variance
    var_P=Qeffectvar/ha2
    var_I=hi2*var_P
    var_E=he2*var_P
    var_e=var_P-Qeffectvar-var_I
    taxa_dupli=NULL
    myY=NULL
    Ieffect_store=NULL
    all_residual=NULL
    ####### simulate Interaction
    ETN.position=sample(1:m,NETN,replace=F)
    SNPE=as.matrix(X[,(ETN.position)])
    trait_name=NULL
    for(t in 1:(NE))
    {
    ####### duplicate the taxa and give different name for each environment    
    taxa_store=paste(taxa,"_",t,sep="")
    taxa_dupli=append(taxa_dupli,taxa_store)
    trait_name=append(trait_name,paste("trait",t,sep=""))
    #interaction
    if(indepent){

    Ieffect=matrix(rnorm(NETN*3,0,1),NETN,3)
    AA_matrix=SNPE
    AT_matrix=SNPE
    TT_matrix=SNPE
    AA_matrix[AA_matrix==1]=2
    AA_matrix[AA_matrix==0]=1
    AA_matrix[AA_matrix==2]=0
    
    AT_matrix[AT_matrix==2]=0
    
    TT_matrix[TT_matrix==1]=0
    TT_matrix[TT_matrix==2]=1
    
    AA=AA_matrix%*%Ieffect[,1]
    AT=AT_matrix%*%Ieffect[,2]
    TT=TT_matrix%*%Ieffect[,3]


    Interaction=AA+AT+TT
    }else{#end indepent
    Ieffect=matrix(rnorm(NETN,0,1),NETN,1)
    Interaction=SNPE%*%Ieffect


    }



    var_I_store=var(Interaction)
    
    bb=as.numeric(var_I/var_I_store)
    Ieffect=sqrt(bb)*Interaction


    Ieffect_store=cbind(Ieffect_store,Ieffect)
    #environment
    Eeffect=rnorm(1,0,sqrt(var_E))
    #residual
    residual1=rnorm(n,0,sqrt(var_e))
    var_e_store=var(residual1)
    cc=as.numeric(var_e/var_e_store)
    residual2=sqrt(cc)*residual1


    #phenotype
    Y[,t]=Qeffect+Ieffect+Eeffect+residual2
    
    #myY=append(myY,Y[,t])
    all_residual=cbind(all_residual,residual2)
    }
    myY=cbind(as.data.frame(taxa),Y)
    colnames(myY)=c("taxa",trait_name)
    
    
    
    
    return(list(Y=myY,u=Qeffect,QTN.position=c(QTN.position,ETN.position),Aeffect=Qeffect,Ieffect=Ieffect_store,e1=all_residual[,1],e2=all_residual[,2]))
} #enf of phenotype simulation function
#=============================================================================================


