rm(list=ls())
memory.limit(102400)
eleMHR<-function(element_mutation_wd,eleMHR_result_wd){
  element_mutation<-read.table(element_mutation_wd,header = T,sep = "\t",stringsAsFactors = F,quote = "")
  id<-unique(element_mutation$gene)
  gene.start<-unique(element_mutation$ElementStart)
  gene.end<-unique(element_mutation$ElementEnd)
  names(gene.start)<-id
  names(gene.end)<-id
  #Divide files into lists by gene
  element_mutation_ls<-split(element_mutation,element_mutation$gene)
  #Run MSEA_cluster
  gene_result<-lapply(element_mutation_ls,MSEA_cluster,gene.start,gene.end,time=1000)
  result<-do.call(rbind.data.frame,gene_result)
  write.table(result,eleMHR_result_wd,sep="\t",quote=FALSE,row.names = F)
}


  ### MSEA_cluster 
  MSEA_cluster= function(gene_data,gene.start,gene.end,time) {
  #####Get a gene start and stop position
  gene.start_1 = gene.start[ gene_data$gene[1] ]
  gene.end_1 =gene.end[ gene_data$gene[1] ]
  ######Get the  mutation location of this gene
  mut_pos = gene_data$start
  
  ######The total number of mutations in this gene
  #######Calculate the length of gene element
  gene.length<-gene.end_1-gene.start_1+1
  
  #####Random
  #####Random 1000 times
  Nm = nrow(gene_data)###The number of mutations in this gene element
  ####Randomly generate mutation sites whose number is the same as the real mutation
  set.seed(123)
  mut_pos.pai = matrix(sample(gene.start_1:gene.end_1, Nm*time, replace=T),nrow=time)
  ### es.random
  es.random<-apply(mut_pos.pai,1,random,gene.length,gene.start_1,gene.end_1)
  
  ### es true
  #####Calculation of real ES
  inc = 1/length(mut_pos)
  ####Total number of mutations per location
  nMut.per.location<-table(gene_data$start)
  dec = 1/(gene.length-length(mut_pos))
  dec.1 = rep(-dec, gene.length)
  names(dec.1)<-as.character(gene.start_1:gene.end_1)
  dec.1[names(nMut.per.location)] = inc * nMut.per.location
  cumsum(dec.1) -> es.cum
  es.true = max(es.cum) - min(es.cum) 
  p = sum(es.random>=es.true)/length(es.random)
  result<-data.frame(gene_data$gene[1],p,nrow(gene_data),names(which.min(es.cum)),min(es.cum),
                     names(which.max(es.cum)),max(es.cum),es.true,stringsAsFactors = F)
  colnames(result)<-c("gene","p","nMut","es_min_pos","es_min","es_max_pos","es_max","es.true")
  return(result)
  } 
  #es.random
  random<-function(mut_data,gene.length,gene.start_1,gene.end_1){
    inc = 1/length(mut_data)
    table(mut_data) -> nMut.per.location
    dec = 1/(gene.length-length(nMut.per.location))
    dec.1 = rep(-dec, gene.length)
    names(dec.1)<-as.character(gene.start_1:gene.end_1)
    dec.1[names(nMut.per.location)] = inc * nMut.per.location
    cumsum(dec.1) -> es.cum
    return(max(es.cum) - min(es.cum))
  }
  
  
 #Run eleMHR
  
  element_mutation_wd<-"/path/data/test_file/element_mutations/CDS.txt"
  eleMHR_result_wd<-"path/data/test_file/R_element_result/CDS_result.txt"
  eleMHR(element_mutation_wd,eleMHR_result_wd)
  
  
  
  
  
  
  
  
  
  
  

