rm(list=ls())
#####For each mutation file, five elements are matched
####Hg19/Hg38 genome element is used according to the reference genome version of the mutation file

mutation_element_match<-function(element_wd,mutation_wd,element_mutation_wd){
  setwd(element_wd) 
  files<-list.files()
  element<-gsub("_hg.*","",files)
  
  #read mutations file
  mutations<-read.table(mutation_wd,header = T,sep = "\t",stringsAsFactors = F,quote = "")
  
  #Divide files into lists by gene
  mutation_ls<-split(mutations,mutations$gene)
  
  
  #####For elements
  i=1
  for (i in 1:length(files)) {
    element_file<-read.table(paste(element_wd,files[i],sep=""),
                             header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
    
    #Delete duplicate gene
    #The minimum and maximum values were taken as the start and end sites of the gene's elements
    gene_start<-aggregate(element_file$START,by=list(element_file$GENE),min)
    gene_end<-aggregate(element_file$END,by=list(element_file$GENE),max)
    element_file<-cbind(gene_start,gene_end$x)
    colnames(element_file)<-c("gene","start","end")
    
    #Match mutations and element
    mutation_file<-mutation_ls$ABCC11
    element_file_fun<-function(mutation_file){
      element_file_index<-which(element_file$gene%in%mutation_file$gene[1])
      if(length(element_file_index)>0){
        element_file_index1<-which(mutation_file$start>=element_file$start[element_file_index] &
                                     mutation_file$end<=element_file$end[element_file_index] )
        if(length(element_file_index1)>0){
          element_file_mut<-mutation_file[element_file_index1,]
          ElementStart<-element_file$start[element_file_index]
          ElementEnd<-element_file$end[element_file_index]
          element_file_mut<-data.frame(ElementStart,ElementEnd,element_file_mut)
          return(element_file_mut)  
        }
        
      }
    }
    element_file_result<-lapply(mutation_ls,element_file_fun)
   
     ####Remove null items
    element_file_result<-element_file_result[vapply(element_file_result, Negate(is.null), NA)]
    
    ##Retain genes with mutations >3
    element_file_nrow<-lapply(element_file_result,nrow)
    element_file_result<-element_file_result[which(element_file_nrow>=4)]
    element_file_result<-do.call(rbind.data.frame,element_file_result)
    
    write.table(element_file_result,paste(element_mutation_wd,element[i],".txt",sep = ""),
                sep="\t",quote=FALSE,row.names = F)
  }
}
###Run mutation_element_match()
element_wd<-"E:/Github/data/element/element_hg38/"
mutation_wd<-"E:/Github/data/test_file/ACC_nonsynonymous_mutation.txt"
element_mutation_wd<-"E:/Github/data/test_file/element_mutations/"
mutation_element_match(element_wd,mutation_wd,element_mutation_wd)




