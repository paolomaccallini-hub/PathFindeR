# file name: PF_Func
#
# Rome June 2025
#
current_dir<-getwd() # current directory
set.seed(12345) # to make results reproducible
#
#-------------------------------------------------------------------------------
# Packages and files
#-------------------------------------------------------------------------------
#
# Load libraries
#
library(dplyr) 
library(rentrez) 
library(httr)
library(jsonlite)
library(curl)
library(biomaRt)
library(stringr)
library(DOSE)
library(stats)
library(rstatix)
library(ReactomePA)
library(igraph)
library(MASS)
library(data.table)
library(clusterProfiler)
library(stats)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(GOplot)
library(readxl) 
library(writexl)
library(ggplot2)
#
#-------------------------------------------------------------------------------
# Build (if necessary) and load STRING database
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/STRING")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Load data for all species, then keep only H. sapiens
#
url<-paste0("https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz")
destfile<-"9606.protein.links.v12.0.txt.gz"
file_path<-file.path(current_dir,"Data/STRING/",destfile)
if(!file.exists(file_path)) {
  print("Downloading STRING database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
url<-paste0("https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz")
destfile<-"9606.protein.info.v12.0.txt.gz"
file_path<-file.path(current_dir,"Data/STRING/",destfile)
if(!file.exists(file_path)) {
  print("Downloading STRING database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
#
print("Loading STRING database...")
#
file.name<-paste0(current_dir,"/Data/STRING/9606.protein.links.v12.0.txt.gz")
STRING.matrix<-read.csv(file.name,sep=" ") # PPI scores
#
file.name<-paste0(current_dir,"/Data/STRING/9606.protein.info.v12.0.txt.gz")
STRING.names<-fread(file.name,sep="\t") # STRING preferred names
colnames(STRING.names)[1]<-"string_protein_id"
#
#-------------------------------------------------------------------------------
# Build (if necessary) and load NCBI database
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/NCBI")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
url<-paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
destfile<-"gene_info.gz"
file_path<-file.path(current_dir,"Data/NCBI/",destfile)
if(!file.exists(file_path)&!file.exists("Data/NCBI/human_gene_info.gz")) {
  test<-0
  while(test==0) {
    print("Downloading NCBI database...")
    attempt<-try(curl_download(url,file_path,quiet=FALSE),silent=T)
    if (class(attempt)!="try-error") {
      test<-1
    }
  }
}   
#
file.name<-"Data/NCBI/human_gene_info.gz"
if (!file.exists(file.name)) {
  print("Editing NCBI database...")
  file.name<-"Data/NCBI/gene_info.gz"
  NCBI.names<-fread(file.name,sep="\t") 
  colnames(NCBI.names)[1]<-"tax_id"
  NCBI.names<-subset.data.frame(NCBI.names,tax_id==9606)
  NCBI.names<-subset.data.frame(NCBI.names,select=c("GeneID","Symbol","Synonyms","description"))
  colnames(NCBI.names)<-c("NCBI.id","name","Synonyms","description")
  file.name<-"Data/NCBI/human_gene_info.gz"
  fwrite(NCBI.names,file.name,sep="\t") 
  file.remove("Data/NCBI/gene_info.gz")
} else {
  print("Loading NCBI database...")
  NCBI.names<-fread(file.name,sep="\t")  
}
gc()
#
#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
#
STRING.co<-0.7 # cut-off for gene interaction in STRING API and STRING data
seed.score<-1 # score for seed genes
#
#-------------------------------------------------------------------------------
# Convert gene symbol to NCBI ID (when possible)
#-------------------------------------------------------------------------------
#
Symbol2NCBI<-function(gene.symbol) {
  #
  # print(paste("Seraching NCBI gene ID for",gene.symbol))
  retries<-100
  for (i in 1:retries) {
    # Use tryCatch to catch errors and warnings
    result<-tryCatch(
      {
        # Perform the search with entrez_search
        NCBI<-lapply(gene.symbol,function(gene) {
          search<-entrez_search(db="gene",term=paste(gene,"[Gene Name] AND human[Organism]"))
        })
      },
      error = function(e) {
        # Check if the error is HTTP 500
        if (grepl("HTTP failure: 500", e$message)) {
          message(paste("Attempt", i, "failed with HTTP 500. Retrying..."))
          return(NA)
        } else {
          stop(e)  # Stop the loop for other errors
        }
      }
    )
    # If result is not NULL, the search was successful
    if (!is.null(result)) {
      if (is.list(result$ids)) {
        return(NA)
      } else if (is.list(NCBI[[1]]$ids[1])) {
        return(NA)
      } else {
        NCBI<-NCBI[[1]]$ids[1]
        return(NCBI)
      }
    }
    # If we've reached this point, the search failed with HTTP 500.
    # Wait before retrying
    Sys.sleep(0.1)
  }
}
#
#-------------------------------------------------------------------------------
# Convert gene symbol to NCBI ID using local data base (not API)
#-------------------------------------------------------------------------------
#
Symbol2NCBI.db<-function(gene.symbol) {
  #
  # Search in "name" first
  #
  df<-subset.data.frame(NCBI.names,name==gene.symbol)
  if (nrow(df)>0) {
    return(as.character(df$NCBI.id[1]))
  } else { # Search among synonyms
    index<-which(grepl(paste0("\\b",gene.symbol,"\\b"),NCBI.names$Synonyms)) 
    if (length(index)==0) {
      return(NA)
    } else {
      return(as.character(NCBI.names$NCBI.id[index[1]]))
    }
  }
}
#
#-------------------------------------------------------------------------------
# Convert gene NCBI ID to gene symbol  (when possible)
#-------------------------------------------------------------------------------
#
NCBI2Symbol<-function(gene.NCBI.id) {
  #
  retries<-100
  for (i in 1:retries) {
    # Use tryCatch to catch errors and warnings
    result<-tryCatch(
      {
        # Perform the search with entrez_search
        entrez_summary(db ="gene",id=gene.NCBI.id)
      },
      error = function(e) {
        # Check if the error is HTTP 500
        if (grepl("HTTP failure: 500", e$message)) {
          message(paste("Attempt", i, "failed with HTTP 500. Retrying..."))
          return(NULL)
        } else {
          stop(e)  # Stop the loop for other errors
        }
      }
    )
    # If result is not NULL, the search was successful
    if (!is.null(result)) {
      if (length(result$name)==0) {
        return(NA)
      } else {
        return(result$name)
      }
    }
    # If we've reached this point, the search failed with HTTP 500.
    # Wait before retrying
    Sys.sleep(0)
  }
}  
#
#-------------------------------------------------------------------------------
# Convert NCBI ID to ENSG ID (when possible)
#-------------------------------------------------------------------------------
#
NCBI2ENSG<-function(NCBI.id) {
  #
  print(paste("Retrieving ENSG ID for", NCBI.id))
  Sys.sleep(1)
  mart <- useEnsembl(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl", 
    mirror = "www",  # Options: "uswest", "useast", "asia", "www"
  )
  result <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
                  filters = "entrezgene_id",
                  values = NCBI.id,
                  mart = mart)
  return(result$ensembl_gene_id[1])
}
#
#-------------------------------------------------------------------------------
# This function asks STRING API its preferred name for a gene
#-------------------------------------------------------------------------------
#
STRING.name<-function(gene) {
  #
  species_id<-9606 # Homo sapiens
  base_url<-"https://string-db.org/api/json/get_string_ids"
  identifiers<-paste(gene,collapse="%0d")
  response<-httr::GET(base_url,query=list(identifiers=identifiers,species=species_id)) 
  #
  if (status_code(response)==200) {
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-jsonlite::fromJSON(response_content) # parse
    if (is.data.frame(data)) {
      preferred.name<-data$preferredName
    } else {
      preferred.name<-NA
    }
  } else {
    stop("Failed to retrieve data. Please check the NCBI IDs and your internet connection.")
  }
  return(preferred.name)
}
#
#-------------------------------------------------------------------------------
# This function find genes that interact with the input gene, 
# according to STRING API
#-------------------------------------------------------------------------------
#
STRING<-function(gene) {
  #
  Sys.sleep(0.1)
  species_id<-9606 # Homo sapiens
  #
  # find STRING's identifier for the inputted gene
  #
  gene_name<-STRING.name(gene)
  if (is.na(gene_name)) gene_name<-"" # necessary, otherwise STRING search gene NA
  #
  # find interacting genes
  #
  base_url<-"https://string-db.org/api/json/network"
  response<-httr::GET(base_url,query=list(identifiers=gene_name,species=species_id,
                                          required_score=STRING.co*1000))
  #
  if (status_code(response)==200) {
    # 
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-jsonlite::fromJSON(response_content) # parse
    #
    if (is.data.frame(data)>0) {
      data<-subset.data.frame(data,score>=STRING.co)
      data<-data.frame(name=data$preferredName_B,score=data$score)
      #
      # we keep for each gene only the entry with highest score
      #
      unique.genes<-unique(data$name) 
      interacting_genes<-data[1:length(unique.genes),]
      for (k in 1:length(unique.genes)) {
        temp.df<-subset.data.frame(data,name==unique.genes[k])
        temp.df<-temp.df[order(temp.df$score,decreasing=T),]
        interacting_genes[k,]<-temp.df[1,]
      }
      return(interacting_genes)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# This function find genes that interact with the input gene, 
# using STRING database (not API)
#-------------------------------------------------------------------------------
#
STRING2<-function(gene) {
  #
  # Find all the interacting genes
  #
  df<-subset.data.frame(STRING.names,preferred_name==gene)
  gene<-df$string_protein_id
  df<-subset.data.frame(STRING.matrix,protein1==gene)
  df<-subset.data.frame(df,combined_score>=STRING.co*1000)
  #
  # edit the output
  #
  df<-df[,-1]
  colnames(df)<-c("name","score")
  #
  if (nrow(df)>0) {
    for (i in 1:nrow(df)) {
      temp.df<-subset.data.frame(STRING.names,string_protein_id==df$name[i])
      df$name[i]<-temp.df$preferred_name
    }
    df$score<-df$score/1000
    return(df)
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# This function find PPI score between two given genes 
# using STRING database (not API)
#-------------------------------------------------------------------------------
#
STRING.db<-function(gene1,gene2) {
  #
  # find PPI score between the two genes
  #
  if (!is.na(gene1)&!is.na(gene2)) {
    # find identifiers
    df<-subset.data.frame(STRING.names,preferred_name==gene1)
    gene1<-df$string_protein_id
    df<-subset.data.frame(STRING.names,preferred_name==gene2)
    gene2<-df$string_protein_id
    # find PPI score
    df<-subset.data.frame(STRING.matrix,protein1==gene1)
    df<-subset.data.frame(df,protein2==gene2)
    if (nrow(df)==1) {
      score<-df$combined_score/1000
    } else {
      score<-0
    }
    if (score<STRING.co) {
      score<-0
    } 
    return(score)
  } else {
    return(NA)
  }
}
# 
#-------------------------------------------------------------------------------
# Expand a gene list using STRING and add NCBI ids
#-------------------------------------------------------------------------------
#
ListExpand<-function(Exper_list) {
  #
  vector<-Exper_list$name
  vector2<-Exper_list$list.name
  vector3<-Exper_list$NCBI.id
  #
  # vector: array of genes to use as seed for PPI
  #
  df<-data.frame(name=vector,score=rep(0,length(vector)), 
                 Predicted=rep("",length(vector)),
                 Predictor=rep("",length(vector)),
                 list.name=rep("",length(vector)))
  #
  for (i in 1:length(vector)) {
    temp<-STRING2(vector[i]) # search for interacting genes
    if (is.data.frame(temp)) {
      temp<-subset.data.frame(temp,name!=vector[i]) # remove seed gene
      temp<-na.omit(temp) # remove NAs
      temp<-cbind(temp,data.frame(Predicted=rep("",nrow(temp))))
      temp<-cbind(temp,data.frame(Predictor=rep(vector[i],nrow(temp))))
      temp$list.name<-rep(vector2[i],nrow(temp))
      df$Predicted[i]<-paste(temp$name,collapse="/") 
      df<-rbind(df,temp) # add to df genes predicted by def$name[i]
    }
  }
  #
  # Add NCBI id and list name
  #
  for (i in 1:nrow(df)) {
    if (i<=length(vector)) {
      df$source[i]<-"seed" # experimental source
      df$score[i]<-seed.score # score for genes of experimental origin
      df$list.name[i]<-vector2[i]
      df$NCBI.id[i]<-vector3[i]
    } else {
      df$source[i]<-"predicted" # predicted from protein-protein interaction
      df$NCBI.id[i]<-Symbol2NCBI.db(df$name[i])
      if (is.na(df$NCBI.id[i])) {
        df$NCBI.id[i]<-Symbol2NCBI(df$name[i]) # try with API
      }
    }
  }
  #
  return(df) # the expanded gene list, including seed genes
}
#
#-------------------------------------------------------------------------------
# For genes that appears more than once, this function generate a single row with
# all the information and with the highest score associated to each gene
#-------------------------------------------------------------------------------
#  
ListCollapse<-function(all.genes.zero) {
  #
  # Indicate predicted genes and predictors
  #
  unique.genes<-unique(all.genes.zero$NCBI.id) # NCBI IDs of unique genes
  all.genes.unique<-all.genes.zero[1:length(unique.genes),]
  for (i in 1:length(unique.genes)) {
    temp<-subset.data.frame(all.genes.zero,NCBI.id==unique.genes[i])
    all.genes.unique$NCBI.id[i]<-temp$NCBI.id[1]
    all.genes.unique$name[i]<-temp$name[1]
    all.genes.unique$source[i]<-paste0(temp$source,collapse="/")
    all.genes.unique$score[i]<-max(as.numeric(temp$score))
    all.genes.unique$list.count[i]<-length(unique(temp$list.name))
    all.genes.unique$list.name[i]<-paste0(unique(temp$list.name),collapse="/")
    all.genes.unique$Predicted[i]<-paste0(temp$Predicted,collapse="/")
    all.genes.unique$Predictor[i]<-paste0(temp$Predictor,collapse="/")
    if (grepl("^(.)(\\1)*$",all.genes.unique$Predicted[i])) {
      all.genes.unique$Predicted[i]<-"" # remove repetitions of "/"
    } 
    if (grepl("^(.)(\\1)*$",all.genes.unique$Predictor[i])) {
      all.genes.unique$Predictor[i]<-"" # remove repetitions of "/"
    }
    temp.str<-all.genes.unique$Predicted[i]
    if (substr(temp.str,1,1)=="/") {
      all.genes.unique$Predicted[i]<-substr(temp.str,2,nchar(temp.str))
    }
    temp.str<-all.genes.unique$Predictor[i]
    if (substr(temp.str,1,1)=="/") {
      all.genes.unique$Predictor[i]<-substr(temp.str,2,nchar(temp.str))
    }
  }
  return(all.genes.unique)
} 
#
#-------------------------------------------------------------------------------
# This function build adjacency matrix from a list of genes using STRING database 
# (not API). Faster than using STRING.db
#-------------------------------------------------------------------------------
#
GeneMatrix<-function(all.genes) {
  #
  print(paste0("I start building gene matrix at time: ",Sys.time()))
  #
  # build a STRING database with only genes of interest and with score above STRING.co
  #
  df1<-STRING.names
  colnames(df1)[2]<-"name"
  df.names<-merge(all.genes,df1,by="name",all.y=F)
  df.names<-subset.data.frame(df.names,select=c("name","string_protein_id"))
  dt.names<-as.data.table(df.names) # data tables should be faster
  remove(df1,df.names)
  #
  df.matrix<-subset.data.frame(STRING.matrix,combined_score>=STRING.co*1000)
  df.matrix<-subset.data.frame(df.matrix,protein1%in%dt.names$string_protein_id)
  df.matrix<-subset.data.frame(df.matrix,protein2%in%dt.names$string_protein_id)
  dt.matrix<-as.data.table(df.matrix) # data tables should be faster
  remove(df.matrix)
  #
  # build gene.matrix
  #
  NG<-nrow(all.genes)
  gene.matrix<-matrix(0,nrow=NG,ncol=NG) 
  rownames(gene.matrix)<-all.genes$name
  colnames(gene.matrix)<-all.genes$name
  #
  time.zero<-as.numeric(Sys.time())
  counter.zero<-0 # pair calculated at time.zero
  todo.counter<-((NG^2)-NG)/2 # total pairs to calculate
  for (i in 1:(NG-1)) {
    #
    gene1<-rownames(gene.matrix)[i]
    #
    # find identifier for gene 1
    #
    dt<-dt.names[name==gene1]
    gene1<-dt$string_protein_id
    #
    # remove data we wont use (dt.matrix is symmetric!)
    #
    dt.matrix<-dt.matrix[protein2!=gene1]
    #
    for (j in (i+1):NG) {
      #
      gene2<-colnames(gene.matrix)[j]
      #
      # find identifier for gene 2
      #
      dt<-dt.names[name==gene2]
      gene2<-dt$string_protein_id
      #
      # find PPI score
      #
      dt<-dt.matrix[protein1==gene1&protein2==gene2]
      #
      if (nrow(dt)==1) {
        score<-dt$combined_score/1000
      } else {
        score<-0
      }
      gene.matrix[i,j]<-score
      #
      # give information on time necessary to finish 
      #
      dice<-sample(seq(1,todo.counter/20,1),1)
      if (dice==50) {
        done.counter<-(j-i)+(2*NG-(i-1+1))*(i-1)/2 # the number of pair calculated
        delta.counter<-done.counter-counter.zero # pairs done from last check
        delta.time<-as.numeric(Sys.time())-time.zero # time from last check (seconds)
        speed<-delta.counter/delta.time # speed for pair calculation
        time.left<-(todo.counter-done.counter)*1/speed # time remaining at current speed (seconds)
        print(paste0("Less than ",round(time.left/60)," minutes to the end, at: ",Sys.time()))
        counter.zero<-done.counter 
        time.zero<-as.numeric(Sys.time())
      }
    }
    #
    # remove data already used
    #
    dt.matrix<-dt.matrix[protein1!=gene1]
    dt.names<-dt.names[string_protein_id!=gene1]
  }
  #
  # add the lower triangle and exit
  #
  gene.matrix[lower.tri(gene.matrix)]<-t(gene.matrix)[lower.tri(gene.matrix)]
  print(paste0("I ended at time: ",Sys.time()))
  return(gene.matrix)
}
#
#-------------------------------------------------------------------------------
# This function select a subset of the adjacency matrix around a specified node.
# A integer number is required to indicate the degree of relationship to be considered.
# Either 1 or 2. It then plots the graph.
#-------------------------------------------------------------------------------
#
SubGeneMatrix<-function(gene.matrix,name,degree) {
  gi<-which(row.names(gene.matrix)==name)
  index<-which(gene.matrix[gi,]!=0)
  if (degree==1) {
    sub.gene.matrix<-gene.matrix[index,index]
  } else if (degree==2) {
    index.1<-index
    for (j in index.1) {
      index<-c(index,which(gene.matrix[j,]!=0))
    }
    index<-unique(index)
    sub.gene.matrix<-gene.matrix[index,index]
  }
  #
  # Build the graph associated with the subset of gene matrix
  #
  graph<-graph_from_adjacency_matrix(sub.gene.matrix,mode="undirected",weighted=TRUE)
  #
  # Color the genes according to the corresponding Expanded Gene List (EGL)
  #
  all.genes<-all.genes[match(rownames(sub.gene.matrix),all.genes$name), ] # correct the order!
  colors<-hcl.colors(nrow(list.matrix),palette="Dark 3",alpha=1)
  vertex_colors<-c()
  for (i in 1:nrow(all.genes)) {
    index<-which(rownames(list.matrix)==all.genes$list.name[i])
    if (length(index)>0) {
      vertex_colors[i]<-colors[index]
    }
  }
  index<-which(is.na(vertex_colors))
  vertex_colors[index]<-"white"
  #
  # Set the color of nodel labels, to highlight seed genes
  #
  node.col<-rep("black",nrow(all.genes))
  for (i in 1:nrow(all.genes)) {
    if("seed"%in%strsplit(all.genes$source[i],"/")[[1]]) {
      node.col[i]<-"red"
    }
  }
  #
  # Plot the image
  #
  tiff(paste0("PF_output/",name,"_",degree,"_graph.tiff"),width=10,height=10,units="in",res=600,compression="lzw")
  plot(graph,
       layout=layout_with_fr(graph,niter=30000,grid="nogrid",dim=2),
       vertex.size=3,
       vertex.label.cex=0.2,
       vertex.frame.color="black",
       vertex.color=vertex_colors,
       vertex.label.color=node.col,
       asp=1,
  )
  legend("topright",                       # or "bottomleft", etc.  
         legend = rownames(list.matrix),   # group names
         col = colors,                     # corresponding colors
         pch = 21,                         # filled circle
         pt.bg = colors,                   # fill color
         pt.cex = 1.5,                     # size of points
         bty = "n")                        # no box around legend
  dev.off()
}
#
#-------------------------------------------------------------------------------
# This function build a matrix where element i-j is the number of genes in common
# between list i and list j. It is the adjacency matrix of a weighted 
# and undirected graph whose nodes are represented by the expanded gene lists
#-------------------------------------------------------------------------------
#
ListMatrix<-function(all.genes) {
  #
  unique.lists<-unique(Exper_list$list.name)
  NL<-length(unique.lists)
  list.matrix<-matrix(data=0,nrow=NL,ncol=NL)
  rownames(list.matrix)<-unique.lists
  colnames(list.matrix)<-unique.lists
  df<-subset.data.frame(all.genes,list.count>1)
  if (nrow(df)>0) {
    for (i in 1:(NL-1)) {
      for (j in (i+1):NL) {
        for (k in 1:nrow(df)) {
          test.vec<-stringr::str_split(df$list.name[k],"/")[[1]]
          if (unique.lists[i]%in%test.vec&unique.lists[j]%in%test.vec) {
            list.matrix[i,j]<-list.matrix[i,j]+1
            list.matrix[j,i]<-list.matrix[j,i]+1
          }
        }
      }  
    }  
  } 
  #
  # Add elements on the diagonal
  #
  if (nrow(all.genes)>0) {
    for (i in 1:NL) {
      for (k in 1:nrow(all.genes)) {
        test.vec<-stringr::str_split(all.genes$list.name[k],"/")[[1]]
        if (unique.lists[i]%in%test.vec) {
          list.matrix[i,i]<-list.matrix[i,i]+1
        }
      }
    }  
  } 
  #
  return(list.matrix)
}
#
#-------------------------------------------------------------------------------
# This function performs Over-representation Analysis over several databases
#-------------------------------------------------------------------------------
#
ORA.fun<-function(all.genes,top.results,list.number) {
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis with KEGG (Kyoto Encyclopedia of Genes and Genomes - 
  # biological pathways)
  #-------------------------------------------------------------------------------
  #
  KEGG_results<-enrichKEGG(gene=all.genes$NCBI.id,organism='hsa',universe=myuniverse)
  KEGG.ORA<-KEGG_results@result
  KEGG.ORA<-KEGG.ORA[,3:ncol(KEGG.ORA)]
  KEGG.ORA$method<-rep("KEGG",nrow(KEGG.ORA))
  #
  matrix.result<-KEGG.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  KEGG.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis with GO (Gene Ontology - functional role of genes within 
  # biological processes, molecular functions, and cellular components)
  #-------------------------------------------------------------------------------
  #
  GO_results<-enrichGO(gene=all.genes$NCBI.id,OrgDb=org.Hs.eg.db,universe=myuniverse)
  GO.ORA<-GO_results@result
  GO.ORA$method<-rep("GO",nrow(GO.ORA))
  #
  matrix.result<-GO.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  GO.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis using Reactome 
  #-------------------------------------------------------------------------------
  #
  pathway_results<-enrichPathway(gene=all.genes$NCBI.id,organism="human",universe=myuniverse)
  Rea.ORA<-pathway_results@result
  Rea.ORA$method<-rep("Reactome",nrow(Rea.ORA))
  #
  matrix.result<-Rea.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  Rea.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis using Disease Ontology (DO)
  #-------------------------------------------------------------------------------
  #
  DO.ORA<-enrichDO(gene=all.genes$NCBI.id,universe=myuniverse)
  DO.ORA<-DO.ORA@result
  DO.ORA$method<-rep("DO",nrow(DO.ORA)) 
  #
  matrix.result<-DO.ORA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  DO.ORA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Merge results from Over-representation Analysis and keep only significant ones
  #-------------------------------------------------------------------------------
  #
  ORA<-rbind(KEGG.ORA,GO.ORA,Rea.ORA,DO.ORA)
  if (nrow(ORA)>0) {
    ORA<-subset.data.frame(ORA,p.adjust<=0.05)
  }
  if (nrow(ORA)>0) {
    ORA<-subset.data.frame(ORA,select=c("ID","Description","GeneRatio",
                                        "pvalue","p.adjust","geneID","method"))
    ORA$gene.name<-rep(NA,nrow(ORA))
    ORA$list.name<-rep(NA,nrow(ORA))
    ORA$list.count<-rep(0,nrow(ORA))
    for (i in 1:nrow(ORA)) {
      vector1<-stringr::str_split(ORA$geneID[i],"/")[[1]]
      vector2<-vector1
      vector3<-vector1
      for (j in 1:length(vector1)) {
        df<-subset.data.frame(all.genes,NCBI.id==as.numeric(vector1[j]))
        vector2[j]<-df$name
        vector3[j]<-df$list.name
      }
      ORA$gene.name[i]<-paste0(unique(vector2),collapse="/")
      ORA$list.name[i]<-paste0(unique(vector3),collapse="/")
      ORA$list.count[i]<-length(vector3)
    }
    for (i in 1:nrow(ORA)) {
      vector1<-unique(stringr::str_split(ORA$list.name[i],"/")[[1]])
      ORA$list.name[i]<-paste0(vector1,collapse="/")
      ORA$list.count[i]<-length(vector1)
    }
  } 
  #
  # Keep only terms with at least two controbuting gene lists
  #
  ORA<-subset.data.frame(ORA,list.count>=list.number)
  #
  # Write gene symbols in alphabetical order
  #
  if (nrow(ORA)>0) {
    for (i in 1:nrow(ORA)) {
      vec<-sort(str_split(ORA$gene.name[i],"/")[[1]])
      vec<-paste0(vec,collapse="/")
      ORA$gene.name[i]<-vec
    }
  } 
  #
  #-------------------------------------------------------------------------------
  # Plots for GO
  #-------------------------------------------------------------------------------
  #
  if (nrow(GO.ORA)>0) {
    index<-which(GO_results@result$ID%in%ORA$ID)
    GO_results@result<-GO_results@result[index,]
    options(ggrepel.max.overlaps = 1000)
    tiff("PF_output/ORA/GO/GO_ORA_goplot.tiff",width=20,height=8,units="in",res=600,compression="lzw")
    p<-goplot(GO_results,showCategory=10)
    print(p)
    dev.off()
    #
    tiff("PF_output/ORA/GO/GO_ORA_cnetplot.tiff",width=20,height=20,units="in",res=600,compression="lzw")
    index<-which(ORA$ID%in%GO_results@result$ID)
    GO_results@result$geneID<-ORA$gene.name[index] # use gene symbols for this plot
    p<-cnetplot(GO_results,showCategory=10,circular=T,colorEdge=T,layout="fr") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggtitle("Top Enriched GO Terms and Genes")
    print(p)
    dev.off()  
  }
  #
  return(ORA)
}
#
#-------------------------------------------------------------------------------
# This function performs gene set enrichment analysis over several databases
#-------------------------------------------------------------------------------
#
GSEA.fun<-function(all.genes,top.results,list.number,score.type) {
  #
  #-------------------------------------------------------------------------------
  # Gene set enrichment analyses (GSEA) on expanded list
  #-------------------------------------------------------------------------------
  #
  # prepare ranked gene list
  #
  geneList<-setNames(as.numeric(all.genes[,score.type]),all.genes$NCBI.id)
  #
  # Identify duplicated ranks (if any)
  #
  duplicated_ranks<-duplicated(geneList)|duplicated(geneList,fromLast=T)
  #
  # Add small random noise only to duplicates
  #
  geneList[duplicated_ranks]<-geneList[duplicated_ranks]+
    runif(sum(duplicated_ranks),-0.0001,0.0001)
  #
  geneList<-sort(geneList,decreasing=TRUE)
  head(geneList)
  #
  #-------------------------------------------------------------------------------
  # KEGG GSEA
  #-------------------------------------------------------------------------------
  #
  KEGG.GSEA<-gseKEGG(geneList=geneList,
                     organism='hsa',
                     scoreType="pos",
                     by="fgsea",
                     seed=T,
                     pvalueCutoff=0.05, 
                     eps=0,
                     verbose=F,
                     nPermSimple=20000)
  KEGG.GSEA<-KEGG.GSEA@result
  #
  matrix.result<-KEGG.GSEA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  KEGG.GSEA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # GO GSEA
  #-------------------------------------------------------------------------------
  #
  GO.GSEA<-gseGO(geneList=geneList,
                 OrgDb=org.Hs.eg.db,
                 scoreType="pos",
                 by="fgsea",
                 seed=T,
                 pvalueCutoff=0.05,
                 eps=0,
                 verbose=F,
                 nPermSimple=20000)
  GO_results<-GO.GSEA # for plotting
  GO.GSEA<-GO.GSEA@result
  #
  matrix.result<-GO.GSEA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  GO.GSEA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Reactome pathway Gene Set Enrichment Analysis
  #-------------------------------------------------------------------------------
  #
  Rea.GSEA<-gsePathway(geneList,
                       organism="human",
                       scoreType="pos",
                       by="fgsea",
                       seed=T,
                       pvalueCutoff=0.05,
                       eps=0,
                       verbose=F,
                       nPermSimple=20000)
  Rea.GSEA<-Rea.GSEA@result
  #
  matrix.result<-Rea.GSEA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  Rea.GSEA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Disease ontology GSEA 
  #-------------------------------------------------------------------------------
  #
  DO.GSEA<-gseDO(geneList,
                 scoreType="pos",
                 ont = "HDO",
                 organism ="hsa",
                 by="fgsea",
                 seed=T,
                 pvalueCutoff=0.05,
                 eps=0,
                 verbose=F,
                 nPermSimple=20000)
  DO.GSEA<-DO.GSEA@result
  #
  matrix.result<-DO.GSEA
  matrix.result<-matrix.result[order(matrix.result$p.adjust),]
  if (nrow(matrix.result)>top.results) matrix.result<-matrix.result[1:top.results,]
  DO.GSEA<-matrix.result
  #
  #-------------------------------------------------------------------------------
  # Merge results from Gene Set Enrichment Analysis and keep only significant ones
  #-------------------------------------------------------------------------------
  #
  GSEA<-rbind(KEGG.GSEA,GO.GSEA,Rea.GSEA,DO.GSEA)
  if (nrow(GSEA)>0) {
    GSEA<-subset.data.frame(GSEA,p.adjust<=0.05)
  }
  if (nrow(GSEA)>0) {
    GSEA<-subset.data.frame(GSEA,select=c("ID","Description","setSize",
                                          "pvalue","p.adjust","core_enrichment"))
    #
    # edit GSEA adding gene names, list names, list count
    #
    GSEA$gene.name<-rep(NA,nrow(GSEA))
    GSEA$list.name<-rep(NA,nrow(GSEA))
    GSEA$list.count<-rep(0,nrow(GSEA))
    for (i in 1:nrow(GSEA)) {
      vector1<-stringr::str_split(GSEA$core_enrichment[i],"/")[[1]]
      vector2<-vector1
      vector3<-vector1
      for (j in 1:length(vector1)) {
        df<-subset.data.frame(all.genes,NCBI.id==as.numeric(vector1[j]))
        vector2[j]<-df$name
        vector3[j]<-unique(df$list.name)
      }
      GSEA$gene.name[i]<-paste0(vector2,collapse="/")
      GSEA$list.name[i]<-paste0(unique(vector3),collapse="/")
    }
    for (i in 1:nrow(GSEA)) {
      vector1<-unique(stringr::str_split(GSEA$list.name[i],"/")[[1]])
      GSEA$list.name[i]<-paste0(vector1,collapse="/")
      GSEA$list.count[i]<-length(vector1)
    }
    #
    # Keep only terms with at least two controbuting gene lists
    #
    GSEA<-subset.data.frame(GSEA,list.count>=list.number)
    #
    # Write gene symbols in alphabetical order
    #
    if (nrow(GSEA)>0) {
      for (i in 1:nrow(GSEA)) {
        vec<-sort(str_split(GSEA$gene.name[i],"/")[[1]])
        vec<-paste0(vec,collapse="/")
        GSEA$gene.name[i]<-vec
      }
    } 
  } 
  #
  #-------------------------------------------------------------------------------
  # Plots for GO
  #-------------------------------------------------------------------------------
  #
  if (nrow(GO.GSEA)) {
    index<-which(GO_results@result$ID%in%GSEA$ID)
    GO_results@result<-GO_results@result[index,]
    options(ggrepel.max.overlaps = 1000)
    tiff("PF_output/GSEA/GO/GO_GSEA_goplot.tiff",width=20,height=8,units="in",res=600,compression="lzw")
    p<-goplot(GO_results,showCategory=10)
    print(p)
    dev.off()
    #
    tiff("PF_output/GSEA/GO/GO_GSEA_cnetplot.tiff",width=20,height=20,units="in",res=600,compression="lzw")
    index<-which(GSEA$ID%in%GO_results@result$ID)
    GO_results@result$geneID<-GSEA$gene.name[index] # use gene symbols for this plot
    p<-cnetplot(GO_results,showCategory=10,circular=T,colorEdge=T,layout="fr") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggtitle("Top Enriched GO Terms and Genes")
    print(p)
    dev.off()  
  }
  #
  return(GSEA)
}
#
#-------------------------------------------------------------------------------
# perform ORA with respect to tissue expression
#-------------------------------------------------------------------------------
# 
Tissue.ORA<-function(all.genes) {
  gs<-GeneSet(geneIds=all.genes$name,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  bk<-GeneSet(geneIds=STRING.names$preferred_name,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  output<-teEnrichment(gs,rnaSeqDataset=1,backgroundGenes=bk) # 1 for HPA, 2 for GTEx
  seEnrichmentOutput<-output[[1]]
  enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  enrichmentOutput$p.ajust<-10^-enrichmentOutput$Log10PValue
  #
  # plot
  #
  tiff("PF_output/ORA/Tissue_ORA.tiff",width=10,height=6,units="in",res=600,compression="lzw")
  p<-ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
                              label = Tissue.Specific.Genes,fill = Tissue))+
    geom_bar(stat = 'identity')+
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
    labs(x='', y = '-LOG10(P-Value)')+
    theme_bw()+
    theme(legend.position='none')+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
            element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          panel.grid.major= element_blank(),panel.grid.minor = element_blank())
  print(p)
  dev.off()
}
#
#-------------------------------------------------------------------------------
# Perform Random walk with restart
#-------------------------------------------------------------------------------
# 
RWR<-function(gene.matrix,Exper_list) {
  #
  # Set back probability
  #
  gamma<-0.7 # back probability
  #
  # Build a data frame with genes
  #
  V<-nrow(gene.matrix) # number of states/genes/nodes
  score.df<-data.frame(name=rownames(gene.matrix))
  #
  # Find seeds 
  #
  score.df$source<-"predicted"
  for (i in 1:V) {
    if (score.df$name[i]%in%Exper_list$name) score.df$source[i]<-"seed"
  }
  #
  # Normalize by column the adjacency matrix and obtain the transition matrix M
  #
  M<-gene.matrix
  for (j in 1:ncol(gene.matrix)) {
    C<-sum(M[,j])
    if (C>0) {
      M[,j]<-M[,j]/C
    } 
  }
  #
  # Set start probability for seed genes, taking weights into account
  #
  score.df$p0<-rep(0,nrow(M)) # starting distribution
  for (i in 1:V) {
    if (score.df$source[i]=="seed") {
      index<-which(Exper_list$name==score.df$name[i])
      score.df$p0[i]<-max(Exper_list$weight[index]) 
    }
  }
  TOT<-sum(score.df$p0)
  score.df$p0<-score.df$p0/TOT
  #
  # Perform random walk with restart using iterative method
  #
  delta<-1
  p1<-score.df$p0
  p2<-p1
  count<-0 # number of iterations
  while (delta>10^-6) {
    p2<-(1-gamma)*p1%*%M+gamma*score.df$p0
    delta<-sum(abs(p2-p1))
    p1<-p2
    count<-count+1
  }
  score.df$p_star_iter<-p2[1,]
  print(paste0("RWR converged in ",count," iterations"))
  #
  # Perform a random walk with analytic solution
  #
  I<-diag(1,V,V) # identity matrix
  score.df$score.RWR<-as.vector(gamma*solve(I-(1-gamma)*M)%*%score.df$p0)
  #
  # Perform random walk with restart using page_rank
  #
  graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
  damping<-1-gamma
  #
  res<-page_rank(
    graph=graph,
    algo="prpack",     
    directed=FALSE,       
    damping=damping,      
    personalized=score.df$p0,       
    weights=E(graph)$weight
  )
  score.df$p_star_pg<-as.vector(res$vector)
  #
  # compare results
  #
  df<-subset.data.frame(score.df,select=c("p_star_iter","score.RWR","p_star_pg"))
  #
  # check for errors 
  #
  er<-c()
  for (i in 1:nrow(df)) {
    max.p<-max(df[i,2:3])
    min.p<-min(df[i,2:3])
    er[i]<-max.p-min.p
  }
  print(paste0("Worst error between analytic solution and Page Rank: ",max(er)))
  print(paste0("Gene with max error: ",score.df$name[which.max(er)]))
  #
  er<-c()
  for (i in 1:nrow(df)) {
    max.p<-max(df[i,c(1,3)])
    min.p<-min(df[i,c(1,3)])
    er[i]<-max.p-min.p
  }
  print(paste0("Worst error between numeric solution and Page Rank: ",max(er)))
  print(paste0("Gene with max error: ",score.df$name[which.max(er)]))
  #
  # Find components
  #
  comps<-components(graph) # connected subgraphs (components)
  #
  # Build a data frame with genes, scores, and components
  #
  for (i in 1:nrow(score.df)) {
    index<-which(score.df$name[i]%in%names(comps$membership))
    score.df$component[i]<-comps$membership[i]
  }
  #
  # Choose the best score 
  #
  score.df$score.RWR<-score.df$p_star_pg
  #
  # Remove columns 
  #
  score.df<-subset.data.frame(score.df,select=-source)
  score.df<-subset.data.frame(score.df,select=-p0)
  score.df<-subset.data.frame(score.df,select=-p_star_iter)
  score.df<-subset.data.frame(score.df,select=-p_star_pg)
  #
  # Return results
  #
  return(score.df)
}
#
#-------------------------------------------------------------------------------
# Count GTEx tissue classes
#-------------------------------------------------------------------------------
#
Tissue_Class_Count<-function(Exper_list) {
  tissue_class <- list(
    PNS = c(
      "Nerve_Tibial"
    ),
    CNS = c(
      "Brain_Anterior_cingulate_cortex_BA24",
      "Brain_Caudate_basal_ganglia",
      "Brain_Cerebellar_Hemisphere",
      "Brain_Cerebellum",
      "Brain_Cortex",
      "Brain_Frontal_Cortex_BA9",
      "Brain_Hippocampus",
      "Brain_Hypothalamus",
      "Brain_Nucleus_accumbens_basal_ganglia",
      "Brain_Putamen_basal_ganglia",
      "Brain_Spinal_cord_cervical_c-1"
    ),
    Cardiovascular = c(
      "Artery_Aorta",
      "Artery_Coronary",
      "Artery_Tibial",
      "Heart_Atrial_Appendage",
      "Heart_Left_Ventricle"
    ),
    Muscle = c(
      "Muscle_Skeletal"
    ),
    Endocrine = c(
      "Adrenal_Gland",
      "Pancreas",
      "Pituitary",
      "Thyroid"
    ),
    Immune = c(
      "Whole_Blood",
      "Cells_EBV-transformed_lymphocytes",
      "Spleen"
    ),
    Digestive = c(
      "Esophagus_Gastroesophageal_Junction",
      "Esophagus_Mucosa",
      "Esophagus_Muscularis",
      "Colon_Sigmoid",
      "Colon_Transverse",
      "Stomach",
      "Small_Intestine_Terminal_Ileum",
      "Liver"
    ),
    Adipose = c(
      "Adipose_Subcutaneous",
      "Adipose_Visceral_Omentum"
    ),
    Respiratory = c(
      "Lung"
    ),
    Reproductive = c(
      "Breast_Mammary_Tissue",
      "Ovary",
      "Prostate",
      "Testis",
      "Uterus",
      "Vagina"
    ),
    Other = c(
      "Minor_Salivary_Gland",
      "Skin_Not_Sun_Exposed_Suprapubic",
      "Skin_Sun_Exposed_Lower_leg",
      "Cells_Cultured_fibroblasts"
    )
  )
  unique_genes<-unique(Exper_list_2$name)
  Exper_list<-data.frame(name=unique_genes,tissue=unique_genes)
  #
  # Collect gene-tissues associations
  #
  for (i in 1:nrow(Exper_list)) {
    df<-subset.data.frame(Exper_list_2,name==unique_genes[i])
    Exper_list$tissue[i]<-paste0(df$Tissues,collapse="/")
  }
  #
  # Remove duplicates
  #
  for (i in 1:nrow(Exper_list)) {
    Exper_list$tissue[i]<-paste0(unique(str_split(Exper_list$tissue[i],pattern="/")[[1]]),collapse="/")
  }
  #
  # Tissue count
  #
  unique_tissues<-paste0(Exper_list$tissue,collapse="/")
  unique_tissues<-unique(str_split(unique_tissues,pattern="/")[[1]])
  v<-unique_tissues
  unique_tissues<-data.frame(tissue=v,tissue_count=v,class=v,class_count=v)
  for (i in 1:nrow(unique_tissues)) {
    index<-grep(unique_tissues$tissue[i],Exper_list$tissue)
    unique_tissues$tissue_count[i]<-length(index)
  }
  unique_tissues$tissue_count<-as.numeric(unique_tissues$tissue_count)
  #
  # Tissue class
  #
  for (i in 1:nrow(unique_tissues)) {
    for (j in 1:length(tissue_class)) {
      if (unique_tissues$tissue[i]%in%tissue_class[[j]]) {
        unique_tissues$class[i]<-names(tissue_class)[j]
      }
    }
  }
  #
  # Tissue class count
  #
  unique_class<-unique(unique_tissues$class)
  for (i in 1:nrow(unique_tissues)) {
    df<-subset.data.frame(unique_tissues,class==unique_tissues$class[i])
    unique_tissues$class_count[i]<-sum(df$tissue_count)
  }
  unique_tissues$class_count<-as.numeric(unique_tissues$class_count)
  unique_tissues<-unique_tissues[order(unique_tissues$tissue_count,decreasing=T),]
  unique_tissues<-unique_tissues[order(unique_tissues$class_count,decreasing=T),]
  #
  # Save results
  #
  write.table(unique_tissues,file="PF_output/Tissues.csv",sep=",")
  #
  return(unique_tissues)
}

