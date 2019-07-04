library(data.table)
library(tidyr)
library(dplyr)
library(Matrix)


mut=fread("somatic_callset_icgc_oscc_26_symec_oscc_14_samples.v1.maf")

#meth.old=fread("for_isi_40pair_swan_nor_mat_EPIC_with_annotation.txt")
#expr.old=fread("fpkm_values_all_genes_from_40_patients_for_ISI.txt")



meth=fread("for_isi_40pair_swan_nor_mat_EPIC_with_annotation_single.txt")
expr=fread("fpkm_values_all_genes_from_40_patients_for_ISI_single_entry_JUNE.txt")


#b=read.csv("Wholegenome_somatic_12patients.csv")
#b1=b[-which(b[,1]=="Unknown"),] # remove sequences that are not annoted genes
#b2=b1[-c(which(b1[,9]=="IGR"),which(b1[,9]=="Intron")),]  # remove mutation in IGR and Intron
#b3_old=b2[-which(b2[,5]=="X"),]  # remove chr X
#meth_old=read.csv("final_methylation_12_patients.csv",sep=",")
#expr_old=read.csv("rnaseq_12patients.csv",sep=",")


b3=mut
fname.1=names(table(b3[,14]))[grep("RADS",names(table(b3[,14])))]
for(i in 1:length(fname.1))
{
    fname.2=NULL
    fname.2= paste("patient",unlist(strsplit(unlist(strsplit(fname.1[i],split="RADS"))[2],split="D")),sep="")
    b3$Tumor_Sample_Barcode=replace(b3$Tumor_Sample_Barcode,which(b3$Tumor_Sample_Barcode==fname.1[i]),fname.2)

}



# MORE THAN 1 ENTRY IN GENE EXPRESSION FOR SOME GENES. WE CONSIDER THE FIRST ENTRY. THE HISTOGRAMS DIFFER SO THEY ARE NOT COMBINED.

# aa=NULL;for(i in 1:length(all_gene)) aa=c(aa,dim(expr[which(expr$gene_short_name==as.character(all_gene[i])),c(seq(44,83,by=1))])[1])



b3.name=NULL;b3.name=unique(unlist(strsplit(as.character(b3[,1]),split="\""))[seq(2,(2*dim(b3)[1]-1),by=2)])

meth.name=NULL;meth.name=unique(unlist(strsplit(as.character(meth[,128]),split="\""))[seq(2,(2*dim(meth)[1]),by=2)])
expr.name=NULL;expr.name=unique(unlist(strsplit(as.character(expr[,1]),split="\""))[seq(2,(2*dim(expr)[1]),by=2)])

all_gene=expr.name[unique(na.omit(match(meth.name[na.omit(match(b3.name,meth.name))],expr.name)))]



p1=p2=p3=NULL
p11=p22=p33=NULL
p111=p222=p333=NULL
h=1
while(h<=length(all_gene))
{
    # FOR NIBMG DATA SENT FOR FIRST TIME
    
    #nor_expr=as.numeric(expr[which(expr$gene_short_name==as.character(all_gene[h]))[1],c(seq(44,83,by=1))])
    #tum_expr=as.numeric(expr[which(expr$gene_short_name==as.character(all_gene[h]))[1],c(seq(4,43,by=1))])

# FOR NIBMG DATA SENT FOR SECOND TIME IN JUNE

    nor_expr=as.numeric(expr[which(expr$gene_short_name==as.character(all_gene[h]))[1],c(seq(43,82,by=1))])
    tum_expr=as.numeric(expr[which(expr$gene_short_name==as.character(all_gene[h]))[1],c(seq(3,42,by=1))])
    
    nor_meth=as.matrix(meth[which(meth$Gene==as.character(all_gene[h])),c(seq(3,81,by=2))])
    tum_meth=as.matrix(meth[which(meth$Gene==as.character(all_gene[h])),c(seq(2,80,by=2))])
    
    
    id=unique(c(which(nor_expr==0),which(tum_expr==0)))
    
    #Y=as.numeric(log(nor_expr)-log(tum_expr))  # changed on 15.06.19
    
    if(length(id)>0 && dim(nor_meth)[1]>1)
    {
        nor_expr=nor_expr[-id]
        tum_expr=tum_expr[-id]
        nor_meth=as.matrix(nor_meth[,-id])
        tum_meth=as.matrix(tum_meth[,-id])
        Y=as.numeric(log(nor_expr)-log(tum_expr))# changed on 15.06.19
        
        #Y=Y[-id] # changed on 19.12.18
    }
    
    if(length(id)>0 && dim(nor_meth)[1]==1)
    {
        nor_expr=nor_expr[-id]
        tum_expr=tum_expr[-id]
        nor_meth=t(as.matrix(nor_meth[-id]))
        tum_meth=t(as.matrix(tum_meth[-id]))
        Y=as.numeric(log(nor_expr)-log(tum_expr))# changed on 15.06.19
        
        #Y=Y[-id] # changed on 19.12.18
    }


    if(length(id)==0)
    {
        nor_expr=nor_expr
        tum_expr=tum_expr
        nor_meth=nor_meth
        tum_meth=tum_meth
        Y=as.numeric(log(nor_expr)-log(tum_expr))# changed on 15.06.19
        
        #Y=Y # changed on 19.12.18
    }
    #  # changed on 19.12.18
    
    if (length(Y)==1)
    {
        p1=c(p1,2); p2=c(p2,2); p3=c(p3,2)
        p11=c(p11,2); p22=c(p22,2); p33=c(p33,2)
        p111=c(p111,2); p222=c(p222,2); p333=c(p333,2)
        cat(h,"\n")
        h=h+1
    }
    if (length(Y)==0)
    {
        p1=c(p1,2); p2=c(p2,2); p3=c(p3,2)
        p11=c(p11,2); p22=c(p22,2); p33=c(p33,2)
        p111=c(p111,2); p222=c(p222,2); p333=c(p333,2)
        cat(h,"\n")
        h=h+1
    }
    
    
    if (length(Y)<=20)
    {
        p1=c(p1,2); p2=c(p2,2); p3=c(p3,2)
        p11=c(p11,2); p22=c(p22,2); p33=c(p33,2)
        p111=c(p111,2); p222=c(p222,2); p333=c(p333,2)
        cat(h,"\n")
        h=h+1
    }
    
    
    if (length(Y)>20)
    {
        b4=b3[which(b3.name==as.character(all_gene[h])),]
        s0=names(table(b3[,14]))
        if(length(id)>0)  s0=s0[-id]
        if(length(id)==0)  s0=s0
        
        
        
        b4.g=b4[which(b4$Reference_Allele=="G"),]
        b4.gc=b4.g[which(b4.g$Tumor_Seq_Allele2 =="C"),]
        b4.ga=b4.g[which(b4.g$Tumor_Seq_Allele2 =="A"),]
        b4.gt=b4.g[which(b4.g$Tumor_Seq_Allele2=="T"),]
        
        
        
        b4.c=b4[which(b4$Reference_Allele=="C"),]
        b4.cg=b4.c[which(b4.c$Tumor_Seq_Allele2 =="G"),]
        b4.ca=b4.c[which(b4.c$Tumor_Seq_Allele2 =="A"),]
        b4.ct=b4.c[which(b4.c$Tumor_Seq_Allele2=="T"),]
        
        
        
        b4.a=b4[which(b4$Reference_Allele=="A"),]
        b4.ag=b4.a[which(b4.a$Tumor_Seq_Allele2 =="G"),]
        b4.ac=b4.a[which(b4.a$Tumor_Seq_Allele2 =="C"),]
        b4.at=b4.a[which(b4.a$Tumor_Seq_Allele2=="T"),]
        
        b4.t=b4[which(b4$Reference_Allele=="T"),]
        b4.tg=b4.t[which(b4.t$Tumor_Seq_Allele2 =="G"),]
        b4.tc=b4.t[which(b4.t$Tumor_Seq_Allele2 =="C"),]
        b4.ta=b4.t[which(b4.t$Tumor_Seq_Allele2=="A"),]
        
        b4.GC_CG=rbind(b4.gc,b4.cg)   # Mutation type: G:C>C:G or C:G>G:C
        b4.GA_CT=rbind(b4.ga,b4.ct)   # Mutation type: G:C>A:T or C:G>T:A
        b4.GT_CA=rbind(b4.gt,b4.ca)   # Mutation type: G:C>T:A or C:G>A:T
        
        b4.AT_TA=rbind(b4.at,b4.ta)   # Mutation type: A:T>T:A or T:A>A:T
        b4.AC_TG=rbind(b4.ac,b4.tg)   # Mutation type: A:T>C:G or T:A>G:C
        b4.AG_TC=rbind(b4.ag,b4.tc)   # Mutation type: A:T>G:C or T:A>C:G
        
        
        
        id0=as.numeric(na.omit(as.numeric(unlist(strsplit(s0,"patient")))))
        id1_GC_CG=which(is.na(match(id0,as.numeric(na.omit(as.numeric(unlist(strsplit(as.character(b4.GC_CG$Tumor_Sample_Barcode),"patient")))))))==FALSE)
        
        GC_CG_mut_in_gene=rep(0,length(s0))
        GC_CG_mut_in_gene[id1_GC_CG]<-1
        
        
        
        id1_GA_CT=which(is.na(match(id0,as.numeric(na.omit(as.numeric(unlist(strsplit(as.character(b4.GA_CT$Tumor_Sample_Barcode),"patient")))))))==FALSE)
        GA_CT_mut_in_gene=rep(0,length(s0))
        GA_CT_mut_in_gene[id1_GA_CT]<-1
        
        id1_GT_CA=which(is.na(match(id0,as.numeric(na.omit(as.numeric(unlist(strsplit(as.character(b4.GT_CA$Tumor_Sample_Barcode),"patient")))))))==FALSE)
        GT_CA_mut_in_gene=rep(0,length(s0))
        GT_CA_mut_in_gene[id1_GT_CA]<-1
        
        id1_AT_TA=which(is.na(match(id0,as.numeric(na.omit(as.numeric(unlist(strsplit(as.character(b4.AT_TA$Tumor_Sample_Barcode),"patient")))))))==FALSE)
        AT_TA_mut_in_gene=rep(0,length(s0))
        AT_TA_mut_in_gene[id1_AT_TA]<-1
        
        id1_AC_TG=which(is.na(match(id0,as.numeric(na.omit(as.numeric(unlist(strsplit(as.character(b4.AC_TG$Tumor_Sample_Barcode),"patient")))))))==FALSE)
        AC_TG_mut_in_gene=rep(0,length(s0))
        AC_TG_mut_in_gene[id1_AC_TG]<-1
        
        id1_AG_TC=which(is.na(match(id0,as.numeric(na.omit(as.numeric(unlist(strsplit(as.character(b4.AG_TC$Tumor_Sample_Barcode),"patient")))))))==FALSE)
        AG_TC_mut_in_gene=rep(0,length(s0))
        AG_TC_mut_in_gene[id1_AG_TC]<-1
        
        A.real=NULL
        if(length(id1_GC_CG)>0)  A.real=cbind(A.real,GC_CG_mut_in_gene)
        if(length(id1_GC_CG)==0)  A.real=A.real
        
        if(length(id1_GA_CT)>0)  A.real=cbind(A.real,GA_CT_mut_in_gene)
        if(length(id1_GA_CT)==0)  A.real=A.real
        
        if(length(id1_GT_CA)>0)  A.real=cbind(A.real,GT_CA_mut_in_gene)
        if(length(id1_GT_CA)==0)  A.real=A.real
        
        if(length(id1_AT_TA)>0)  A.real=cbind(A.real,AT_TA_mut_in_gene)
        if(length(id1_AT_TA)==0)  A.real=A.real
        
        if(length(id1_AC_TG)>0)  A.real=cbind(A.real,AC_TG_mut_in_gene)
        if(length(id1_AC_TG)==0)  A.real=A.real
        
        if(length(id1_AG_TC)>0)  A.real=cbind(A.real,AG_TC_mut_in_gene)
        if(length(id1_AG_TC)==0)  A.real=A.real
        
        if(dim(nor_meth)[1]>1)
        {
            z1=as.numeric(apply(nor_meth,2,mean))-as.numeric(apply(tum_meth,2,mean))  # uses mean
            z2=as.numeric(apply(nor_meth,2,median))-as.numeric(apply(tum_meth,2,median))  # uses median
            
            z3=as.numeric(apply(nor_meth,2,mean)/apply(nor_meth,2,sd))-as.numeric(apply(tum_meth,2,mean)/apply(tum_meth,2,sd)) # uses standardised mean
            
        }
        
        if(dim(nor_meth)[1]==1)     z1=z2=as.numeric(nor_meth-tum_meth)
    
        A.real=cbind(rep(1,length(s0)),A.real)
        if(dim(as.matrix(A.real))[2]==1)
        
        {
            p1=c(p1,2)
            p2=c(p2,2)
            p3=c(p3,2)
            p11=c(p11,2)
            p22=c(p22,2)
            p33=c(p33,2)
            p111=c(p111,2)
            p222=c(p222,2)
            p333=c(p333,2)
            cat(h,"\n")
            h=h+1
            
        }
        if(dim(A.real)[2]>1)
        {
            s=svd(t(A.real)%*%A.real)
            SV=s$v
            SU=s$u
            
            diag_element=s$d
            id_diag=which(diag_element==0)
            
            if (length(id_diag)>0)
            {
                diag_element=diag_element[-id_diag]
                SV=SV[-id_diag,-id_diag]
                SU=SU[-id_diag,-id_diag]
                A.real=A.real[,-id_diag]
            }
            
            if (length(id_diag)==0)
            {
                diag_element=diag_element
                SV=SV
                SU=SU
                A.real=A.real
            }
            
            if((length(s0)-as.numeric(rankMatrix(as.matrix(A.real)))-1)>0)
            {
                D <- diag(diag_element)
                DD<- 1/(diag_element)
                
                
                AA.real_inv=(SV)%*%diag(DD)%*%t(SU)
                P.real=diag(1,length(s0))- (A.real%*% AA.real_inv %*%t(A.real))
                #For z1.sim (mean of methylation values)
                B1.real=(as.matrix(z1)%*%t(as.matrix(z1)))/as.numeric((t(as.matrix(z1)) %*% P.real %*% as.matrix(z1)))
                V1.real=P.real%*%(diag(1,length(s0))-(B1.real%*%P.real))
                #For z2.sim (median of methylation values)
                B2.real=(as.matrix(z2)%*%t(as.matrix(z2)))/as.numeric((t(as.matrix(z2)) %*% P.real %*% as.matrix(z2)))
                V2.real=P.real%*%(diag(1,length(s0))-(B2.real%*%P.real))
                #For z3.sim (standardised mean of methylation values)
                if(dim(nor_meth)[1]>1)
                {
                    B3.real=(as.matrix(z3)%*%t(as.matrix(z3)))/as.numeric((t(as.matrix(z3)) %*% P.real %*% as.matrix(z3)))
                    V3.real=P.real%*%(diag(1,length(s0))-(B3.real%*%P.real))
                }
                
                
                F_realdata_H01_mean=((sum((Y-mean(Y))^2)-(t(Y)%*%V1.real%*%Y))/(as.numeric(rankMatrix(A.real))))/((t(Y)%*%V1.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                F_realdata_H01_median=((sum((Y-mean(Y))^2)-(t(Y)%*%V2.real%*%Y))/(as.numeric(rankMatrix(A.real))))/((t(Y)%*%V2.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                
                p1=c(p1,pf(F_realdata_H01_mean,as.numeric(rankMatrix(A.real)),length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                p2=c(p2,pf(F_realdata_H01_median,as.numeric(rankMatrix(A.real)),length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                
                if(dim(nor_meth)[1]>1)
                {
                    F_realdata_H01_stdmean=((sum((Y-mean(Y))^2)-(t(Y)%*%V3.real%*%Y))/(as.numeric(rankMatrix(A.real))))/((t(Y)%*%V3.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                    p3=c(p3,pf(F_realdata_H01_stdmean,as.numeric(rankMatrix(A.real)),length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                    
                }
                if(dim(nor_meth)[1]==1) p3=c(p3,2)
                
                beta.cap_realdata.z1=sum((z1-mean(z1))*(Y-mean(Y)))/sum((z1-mean(z1))^2)
                beta.cap_realdata.z2=sum((z2-mean(z2))*(Y-mean(Y)))/sum((z2-mean(z2))^2)
                
                if(dim(nor_meth)[1]>1)
                beta.cap_realdata.z3=sum((z3-mean(z3))*(Y-mean(Y)))/sum((z3-mean(z3))^2)
                
                
                
                F_realdata_H02_mean=(((sum((Y-mean(Y))^2)-((beta.cap_realdata.z1)^2)*sum((z1-mean(z1))^2)-t(Y)%*%V1.real%*%Y))/(as.numeric(rankMatrix(A.real))-1))/((t(Y)%*%V1.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                F_realdata_H02_median=(((sum((Y-mean(Y))^2)-((beta.cap_realdata.z2)^2)*sum((z2-mean(z2))^2)-t(Y)%*%V2.real%*%Y))/(as.numeric(rankMatrix(A.real))-1))/((t(Y)%*%V2.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                
                p11=c(p11,pf(F_realdata_H02_mean,as.numeric(rankMatrix(A.real))-1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                p22=c(p22,pf(F_realdata_H02_median,as.numeric(rankMatrix(A.real))-1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                
                if(dim(nor_meth)[1]>1)
                {
                    F_realdata_H02_stdmean=(((sum((Y-mean(Y))^2)-((beta.cap_realdata.z3)^2)*sum((z3-mean(z3))^2)-t(Y)%*%V3.real%*%Y))/(as.numeric(rankMatrix(A.real))-1))/((t(Y)%*%V1.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                    p33=c(p33,pf(F_realdata_H02_stdmean,as.numeric(rankMatrix(A.real))-1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                }
                
                if(dim(nor_meth)[1]==1) p33=c(p33,2)
                
                
                
                #F_realdata_H02_mean=((sum((Y-mean(Y))^2)-mean((z1-mean(z1))*(Y-mean(Y))))-(t(Y)%*%V1.real%*%Y))/(t(Y)%*%V1.real%*%Y)
                #F_realdata_H02_median=((sum((Y-mean(Y))^2)-mean((z2-mean(z2))*(Y-mean(Y))))-(t(Y)%*%V2.real%*%Y))/(t(Y)%*%V2.real%*%Y)
                #F_realdata_H02_stdmean=((sum((Y-mean(Y))^2)-mean((z3-mean(z3))*(Y-mean(Y))))-(t(Y)%*%V3.real%*%Y))/(t(Y)%*%V3.real%*%Y)
                
                
                
                
                #V1.real.H03=t(P.real+(A.real%*% AA.real_inv %*%t(A.real))%*%B1.real%*%P.real)%*%(P.real+(A.real%*% AA.real_inv %*%t(A.real))%*%B1.real%*%P.real)
                #V2.real.H03=t(P.real+(A.real%*% AA.real_inv %*%t(A.real))%*%B2.real%*%P.real)%*%(P.real+(A.real%*% AA.real_inv %*%t(A.real))%*%B2.real%*%P.real)
                #V3.real.H03=t(P.real+(A.real%*% AA.real_inv %*%t(A.real))%*%B3.real%*%P.real)%*%(P.real+(A.real%*% AA.real_inv %*%t(A.real))%*%B3.real%*%P.real)
                
                F_realdata_H03_mean=((t(Y)%*%P.real%*%Y)-(t(Y)%*%V1.real%*%Y))/((t(Y)%*%V1.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                
                F_realdata_H03_median=((t(Y)%*%P.real%*%Y)-(t(Y)%*%V2.real%*%Y))/((t(Y)%*%V2.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                p111=c(p111,pf(F_realdata_H03_mean,1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                
                p222=c(p222,pf(F_realdata_H03_median,1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                
                if(dim(nor_meth)[1]>1)
                {
                    F_realdata_H03_stdmean=((t(Y)%*%P.real%*%Y)-(t(Y)%*%V3.real%*%Y))/((t(Y)%*%V3.real%*%Y)/(length(s0)-as.numeric(rankMatrix(A.real))-1))
                    p333=c(p333,pf(F_realdata_H03_stdmean,1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE))
                }
                if(dim(nor_meth)[1]==1) p333=c(p333,2)
                
                #cat(as.numeric(pf(F_realdata_H01_mean,as.numeric(rankMatrix(A.real)),length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t \t")
                #cat(as.numeric(pf(F_realdata_H01_median,as.numeric(rankMatrix(A.real)),length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H01_stdmean,as.numeric(rankMatrix(A.real)),length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H02_mean,as.numeric(rankMatrix(A.real))-1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H02_median,as.numeric(rankMatrix(A.real))-1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H02_stdmean,as.numeric(rankMatrix(A.real))-1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H03_mean,1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H03_median,1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                #cat(as.numeric(pf(F_realdata_H03_stdmean,1,length(s0)-as.numeric(rankMatrix(A.real))-1,lower.tail = FALSE)),"\t\t")
                
                cat(h,"\n")
                h=h+1
            }
            
            if((length(s0)-as.numeric(rankMatrix(as.matrix(A.real)))-1)==0)
            {
                p1=c(p1,2)
                p2=c(p2,2)
                p3=c(p3,2)
                p11=c(p11,2)
                p22=c(p22,2)
                p33=c(p33,2)
                p111=c(p111,2)
                p222=c(p222,2)
                p333=c(p333,2)
                cat(h,"\n")
                h=h+1
            }
            if((length(s0)-as.numeric(rankMatrix(as.matrix(A.real)))-1)<0)
            {
                p1=c(p1,2)
                p2=c(p2,2)
                p3=c(p3,2)
                p11=c(p11,2)
                p22=c(p22,2)
                p33=c(p33,2)
                p111=c(p111,2)
                p222=c(p222,2)
                p333=c(p333,2)
                cat(h,"\n")
                h=h+1
            }
            
        }
        
        
    }
    
}
result.1.2=cbind(p1,p11,p111,p2,p22,p222)
colnames(result.1.2)=paste(c("p.H01.mean","p.H02.mean","p.H03.mean","p.H01.median","p.H02.median","p.H03.median"))

result.3=cbind(p3,p33,p333)
colnames(result.3)=paste(c("p.H01.stdmean","p.H02.stdmean","p.H03.stdmean"))

fname.1.2=paste("pvalues_unadjusted_mean_median.txt")
fname.3=paste("pvalues_unadjusted_stdmean.txt")

write.table(result.1.2,fname.1.2,row.names=T,col.names=T,sep="\t",quote=F)
write.table(result.3,fname.3,row.names=T,col.names=T,sep="\t",quote=F)


#p1.2_tab=read.table("p1.2.txt",sep="\t")
#p1=p1.2_tab[,1];p11=p1.2_tab[,2];p111=p1.2_tab[,3]
#p2=p1.2_tab[,4];p22=p1.2_tab[,5];p222=p1.2_tab[,6]

ide=1:length(all_gene)
ide1=ide[-which(p1==2)]
temp=matrix(2,ncol=6,nrow=length(all_gene))
p1.adj=p.adjust(p1[-which(p1==2)],"BH")
p11.adj=p.adjust(p11[-which(p11==2)],"BH")
p111.adj=p.adjust(p111[-which(p111==2)],"BH")
p2.adj=p.adjust(p2[-which(p2==2)],"BH")
p22.adj=p.adjust(p22[-which(p22==2)],"BH")
p222.adj=p.adjust(p222[-which(p222==2)],"BH")

temp[ide1,]=cbind(p1.adj,p11.adj,p111.adj,p2.adj,p22.adj,p222.adj)
rownames(temp)[ide1]=all_gene[ide1]

#p3_tab=read.table("p3.txt",sep="\t")
#p3=p3_tab[,1];p33=p3_tab[,2];p333=p3_tab[,3]

ide2=ide[-which(p3==2)]
temp1=matrix(2,ncol=3,nrow=length(all_gene))
p3.adj=p.adjust(p3[-which(p3==2)],"BH")
p33.adj=p.adjust(p33[-which(p33==2)],"BH")
p333.adj=p.adjust(p333[-which(p333==2)],"BH")

temp1[ide2,]=cbind(p3.adj,p33.adj,p333.adj)
rownames(temp1)[ide2]=all_gene[ide2]


fname.adj.1.2=paste("pvalues_adjusted_mean_median.txt")
fname.adj.3=paste("pvalues_adjusted_stdmean.txt")


result.adj.1.2=temp[ide1,]
result.adj.3=temp1[ide2,]

colnames(result.adj.1.2)=paste(c("p.H01.mean","p.H02.mean","p.H03.mean","p.H01.median","p.H02.median","p.H03.median"))

colnames(result.adj.3)=paste(c("p.H01.stdmean","p.H02.stdmean","p.H03.stdmean"))


write.table(result.adj.1.2,fname.adj.1.2,row.names=T,col.names=T,sep="\t",quote=F)
write.table(result.adj.3,fname.adj.3,row.names=T,col.names=T,sep="\t",quote=F)



##################################
########                  ########
########     RESULTS      ########
########                  ########
##################################

#####  as.character(all_gene[which(temp1[,1]<0.05)])  # CECR1 gene: H01 is rejected using standardised mean
##### as.character(all_gene[which(temp1[,3]<0.05)])   # CECR1 gene: H03 is rejected using standardised mean




















