# file name: PathFindeR_Modules
#
# Latest update: 22nd September 2025
#
source("PathFindeR_Func.R")
#
#-----------------------------------------------------------------------------
# Read and edit input gene list (seeds)
#-----------------------------------------------------------------------------
#
# Read 32 candidate genes form DecodeME preprint
#
Exper_list_1<-read.csv("My_genes_DecodeME_preprint.csv",sep=";",header=T) 
#
# Read candidate genes from fine-mapping
#
Exper_list_2<-read.csv("My_genes_DecodeME.csv",sep=";",header=T) # experimental data (UK Biobank)
#
# Subset the 32 candidate genes
#
index<-which(Exper_list_1$name%in%Exper_list_2$name)
Exper_list<-Exper_list_1[index,]
Exper_list<-subset.data.frame(Exper_list,select=c("name","NCBI.id","list.name","weight","risk.locus","region"))
#
# Check NCBI id
#
for (i in 1:nrow(Exper_list)) {
  if (Exper_list$NCBI.id[i]!=Symbol2NCBI.db(Exper_list$name[i])) {
    print(paste("There is a mismatch with NCBI id of:",Exper_list$name[i]))
    Exper_list$NCBI.id[i]<-Symbol2NCBI.db(Exper_list$name[i])
  }
}
#
# replace gene name with STRING's preferred name
#
gene.wo.STRING<-0
for (i in 1:nrow(Exper_list)) {
  new.name<-STRING.name(Exper_list$name[i])
  if (!is.na(new.name)) {
    Exper_list$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Exper_list$name[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
  }
}
print(paste(gene.wo.STRING,"genes are not present in STRING's database"))
#
#-------------------------------------------------------------------------------
# Create output folder, if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/ORA")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/ORA/KEGG")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/ORA/Reactome")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/ORA/GO")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/GSEA")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/GSEA/KEGG")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/GSEA/Reactome")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
current_dir<-getwd()
folder_path<-file.path(current_dir,"PF_output/GSEA/GO")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
#
#-------------------------------------------------------------------------------
# Prepare background genes
#-------------------------------------------------------------------------------
#
file.name<-"Data/My_Universe.tsv"
if (!file.exists(file.name)) {
  print("Building backgroud for Over-representation analysis. This may take a while...")
  STRING.names$NCBI.id<-rep(NA,nrow(STRING.names))
  for (i in 1:nrow(STRING.names)) {
    STRING.names$NCBI.id[i]<-Symbol2NCBI.db(STRING.names$preferred_name[i])
  }
  myuniverse<-as.character(STRING.names$NCBI.id)  
  write.table(myuniverse,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
} else {
  myuniverse<-fread(file.name)
  myuniverse<-as.character(myuniverse$x)
}
#
#-----------------------------------------------------------------------------
# Expand list with STRING 
#-----------------------------------------------------------------------------
#
all.genes.zero<-ListExpand(Exper_list)
all.genes.zero<-na.omit(all.genes.zero)
#
#-------------------------------------------------------------------------------
# Build a list with a single row for each gene
#-------------------------------------------------------------------------------
#
all.genes<-ListCollapse(all.genes.zero)
#
#-------------------------------------------------------------------------------
# Build the adjacency matrix associated with the merged gene list
#-------------------------------------------------------------------------------
#
gene.matrix<-GeneMatrix(all.genes)
#
#-------------------------------------------------------------------------------
# Use random walk with restart for scoring and component assignation
#-------------------------------------------------------------------------------
#
score.df<-RWR(gene.matrix,Exper_list) # standard RWR with also component analysis
all.genes<-merge(all.genes,score.df,by="name",all.x=T)
#
# Add rank with respect to RWR.raw 
#
all.genes<-all.genes[order(all.genes$score.RWR,decreasing=T),]
all.genes$rank.raw<-seq(1,nrow(all.genes),1) # add a rank
#
# Add cumulative stationary probability (CSP)  with respect to RWR.raw
#
all.genes<-all.genes[order(all.genes$score.RWR,decreasing=T),]
all.genes$CSP<-cumsum(all.genes$score.RWR)
#
# Save a copy of the results
#
file.name<-"PF_output/All_genes.tsv"
write.table(all.genes,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# Build a graph with the expanded gene lists 
#-------------------------------------------------------------------------------
#
list.matrix<-ListMatrix(all.genes)
#
#-------------------------------------------------------------------------------
# Save relevant results as RSD file
#-------------------------------------------------------------------------------
#
Disease_results_list<-list(Exper_list,all.genes,gene.matrix,list.matrix)
saveRDS(Disease_results_list,file="PF_output/Disease_results_list.rds")
#
#-------------------------------------------------------------------------------
# Read relevant results as RSD file
#-------------------------------------------------------------------------------
#
Disease_results_list<-readRDS(file="PF_output/Disease_results_list.rds")
Exper_list<-Disease_results_list[[1]]
all.genes<-Disease_results_list[[2]]
gene.matrix<-Disease_results_list[[3]]
list.matrix<-Disease_results_list[[4]]
#
#-------------------------------------------------------------------------------
# Build the graph associated with the Merged Gene List (MGL) and plot it
#-------------------------------------------------------------------------------
#
graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
#
# Color the genes according to the corresponding Expanded Gene List (EGL)
#
all.genes<-all.genes[match(rownames(gene.matrix),all.genes$name), ] # correct the order!
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
tiff("PF_output/All_genes_graph.tiff",width=10,height=10,units="in",res=600,compression="lzw")
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
#
# Save graph for cytoscape
#
edges<-as.data.frame(as_edgelist(graph))
weights<-E(graph)$weight
edges_df<-data.frame(source=edges[,1],target=edges[,2],interaction="interacts_with",
                     weight=weights)
file_name="PF_output/All_genes_cytoscape.tsv"
write.table(edges_df,file=file_name,sep="\t",row.names=FALSE,quote=FALSE)
#
#-------------------------------------------------------------------------------
# Plot the graph associated with the expanded gene lists (EGLs)
#-------------------------------------------------------------------------------
#
# For each overlap different from zero, calculate a p value by hypergeometric test
#
NL<-nrow(list.matrix)
NB<-nrow(STRING.names) # total background (white and black balls: m+n)
list.matrix.pvalue<-list.matrix
for (i in 1:(NL-1)) {
  for (j in (i+1):NL) {
    list.matrix[i,j]<-as.numeric(list.matrix[i,j])
    if (list.matrix[i,j]==0) {
      list.matrix.pvalue[i,j]<-1
    } else {
      N1<-list.matrix[i,i] # size of EPF_i (number of white balls, m)
      N2<-list.matrix[j,j] # size of EPF_j (number of extracted balls, k)
      NR<-list.matrix[i,j] # size of overlap (number of white balls extracted, q) 
      list.matrix.pvalue[i,j]<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F)) # P(X>NR-1)=P(X>=NR)
      list.matrix.pvalue[i,j]<-as.numeric(formatC(list.matrix.pvalue[i,j],format="e",digits=2))
    }
    list.matrix.pvalue[j,i]<-list.matrix.pvalue[i,j]
  }
}
for (i in 1:NL) {
  list.matrix.pvalue[i,i]<-NA
}
#
# Generate graph from list.matrix
#
temp.matrix<-list.matrix
diag(temp.matrix)<-0 # this is necessary to avoid loops
graph<-graph_from_adjacency_matrix(temp.matrix,mode="undirected",weighted=T)
#
# Generate graph from list.matrix.pvalue
#
temp.matrix<-list.matrix.pvalue
diag(temp.matrix)<-0 # this is necessary to avoid loops
graph.pvalue<-graph_from_adjacency_matrix(temp.matrix,mode="undirected",weighted=T)
#
# Prepare colors for the EGLs
#
colors<-hcl.colors(nrow(list.matrix),palette="Dark 3",alpha=1)
#
# Prepare width for edges
#
edge.w<-E(graph)$weight
edge.w.scaled<-(edge.w-min(edge.w))/(max(edge.w)-min(edge.w))*8+5
#
# Plot it
#
tiff("PF_output/All_genes_list_matrix.tiff",width=10,height=10,units="in",res=600,compression="lzw")
plot(graph,
     vertex.size=30,
     vertex.label.cex=2.,
     vertex.label=paste0(rownames(list.matrix),"\n(size = ",diag(list.matrix),")"),
     vertex.frame.color="black",
     vertex.color=colors,
     edge.label=paste0("overlap: ",E(graph)$weight,"\n p = ",E(graph.pvalue)$weight),
     edge.label.cex=1.5,
     edge.width=edge.w.scaled,
     asp=1,
     )
dev.off()
#
#-------------------------------------------------------------------------------
# Prepare and load disease modules from relevant literature
#-------------------------------------------------------------------------------
#
# Read ME/CFS module form Zhang S. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)
# We use supplementary table 2.
#
Zhang_module<-read_xlsx("media-9.xlsx") 
Zhang_module<-subset.data.frame(Zhang_module,q_value<0.02) # 115 genes associated with ME/CFS  
#
# Use STRING preferred name for genes
#
gene.wo.STRING<-0
index<-c()
Zhang_module$name<-rep(NA,nrow(Zhang_module))
for (i in 1:nrow(Zhang_module)) {
  new.name<-STRING.name(Zhang_module$Gene[i])
  if (!is.na(new.name)) {
    Zhang_module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Zhang_module$Gene[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) Zhang_module<-Zhang_module[-index,] # only genes in STRING
#
# Read and edit Open Target Platforms associations for ME/CFS 
#
OT_file<-list.files(pattern="-associated-targets-",full.names=T)
OT_module<-fread(OT_file) 
#
# Use STRING preferred name
#
index<-c()
gene.wo.STRING<-0
for (i in 1:nrow(OT_module)) {
  new.name<-STRING.name(OT_module$symbol[i])
  if (!is.na(new.name)) {
    OT_module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(OT_module$symbol[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) OT_module<-OT_module[-index,]
#
# Comparison with proteomic study (https://pmc.ncbi.nlm.nih.gov/articles/PMC12254397)
# Download: https://pmc.ncbi.nlm.nih.gov/articles/instance/12254397/bin/44321_2025_258_MOESM7_ESM.xlsx
# Put the file in the main folder.
#
# Edit file: keep only STRING genes
#
Results<-list()
if (!file.exists("Data/UKBB_proteomics.xlsx")) {
  for (s in c(2,3)) {
    #
    Proteo_module<-read_xlsx("44321_2025_258_MOESM7_ESM.xlsx",sheet=s) 
    #
    # Use STRING preferred name
    #
    index<-c()
    gene.wo.STRING<-0
    for (i in 1:nrow(Proteo_module)) {
      new.name<-STRING.name(Proteo_module$name[i])
      if (!is.na(new.name)) {
        Proteo_module$name[i]<-new.name
      } else if (is.na(new.name)) {
        print(paste(Proteo_module$name[i],"is not present in STRING's database"))
        gene.wo.STRING<-gene.wo.STRING+1
        index<-c(index,i)
      }
    }
    if (length(index)>0) {
      Results[[s-1]]<-Proteo_module[-index,] # remove genes not in STRING
    } else {
      Results[[s-1]]<-Proteo_module
    }
  }
  write_xlsx(Results,"Data/UKBB_proteomics.xlsx",col_names=T)   
}
#
#-------------------------------------------------------------------------------
# Prepare different candidate disease modules
#-------------------------------------------------------------------------------
#
Disease.module_list<-list()
Disease.module_list[[1]]<-all.genes
names(Disease.module_list)[1]<-"All genes"
k<-2
for (CSPco in c(0.80,0.85,0.90,0.95)) {
  all.genes<-all.genes[order(all.genes$score.RWR,decreasing=T),]
  n.top<-which(all.genes$CSP>=CSPco)[1]
  Disease.module_list[[k]]<-all.genes[1:n.top,]
  names(Disease.module_list)[k]<-paste0("Disease.module.raw.",CSPco)
  k<-k+1
}
#
#-------------------------------------------------------------------------------
# Hypergeometric tests for comparison with relevant gene modules 
#-------------------------------------------------------------------------------
#
# Prepare a list for the results
#
myhypertest_list<-list()
#
# Set A and test counter
#
test_count<-0
mygenes<-Disease.module_list
#
# Hypergeometric test against Zhang module
#
setB<-Zhang_module
setB.name<-"Zhang S 2025"
NB<-15000
myhypertest<-data.frame(setA=names(mygenes),setB=rep(setB.name,length(mygenes)))
for (i in 1:length(mygenes)) {
  intersection<-merge(mygenes[[i]],setB,by="name")  
  N1<-nrow(mygenes[[i]])
  N2<-nrow(setB)
  NR<-nrow(intersection)
  pvalue<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F))
  myhypertest$sizeA[i]<-N1
  myhypertest$sizeB[i]<-N2
  myhypertest$overlap[i]<-NR
  myhypertest$background[i]<-NB
  myhypertest$pvalue[i]<-formatC(pvalue,format="e",digits=2)
  myhypertest$genes[i]<-paste0(intersection$name,collapse="/")
}
test_count<-test_count+1
myhypertest_list[[test_count]]<-myhypertest
#
# Hypergeometric test against Open Targets
#
NB<-15000
setB<-OT_module
setB.name<-"Open Targets"
myhypertest<-data.frame(setA=names(mygenes),setB=rep(setB.name,length(mygenes)))
for (i in 1:length(mygenes)) {
  intersection<-merge(mygenes[[i]],setB,by="name")  
  N1<-nrow(mygenes[[i]])
  N2<-nrow(setB)
  NR<-nrow(intersection)
  pvalue<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F))
  myhypertest$sizeA[i]<-N1
  myhypertest$sizeB[i]<-N2
  myhypertest$overlap[i]<-NR
  myhypertest$background[i]<-NB
  myhypertest$pvalue[i]<-formatC(pvalue,format="e",digits=2)
  myhypertest$genes[i]<-paste0(intersection$name,collapse="/")
}
test_count<-test_count+1
myhypertest_list[[test_count]]<-myhypertest
#
# Hypergeometric test against proteomics
#
universe<-read_xlsx("Data/UKBB_proteomics.xlsx",sheet=1) # total effect 
NB<-nrow(universe)
setB<-subset.data.frame(universe,BH_pval_combined_te<=0.05)  
setB.name<-"Proteomics"
myhypertest<-data.frame(setA=names(mygenes),setB=rep(setB.name,length(mygenes)))
for (i in 1:length(mygenes)) {
  setA<-merge(mygenes[[i]],universe,by="name") # keep only genes in universe
  intersection<-merge(setA,setB,by="name")  
  N1<-nrow(setA)
  N2<-nrow(setB)
  NR<-nrow(intersection)
  pvalue<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F))
  myhypertest$sizeA[i]<-N1
  myhypertest$sizeB[i]<-N2
  myhypertest$overlap[i]<-NR
  myhypertest$background[i]<-NB
  myhypertest$pvalue[i]<-formatC(pvalue,format="e",digits=2)
  myhypertest$genes[i]<-paste0(intersection$name,collapse="/")
}
test_count<-test_count+1
myhypertest_list[[test_count]]<-myhypertest
#
# Build the final data frame and save it
#
df<-myhypertest_list[[1]]
for (i in 2:length(myhypertest_list)) {
  df<-rbind(df,myhypertest_list[[i]])
}
#
file.name<-"PF_output/replication.tsv"
write.table(df,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# Select a disease module 
#-------------------------------------------------------------------------------
#
all.genes<-all.genes[order(all.genes$score.RWR,decreasing=T),]
n.top<-which(all.genes$CSP>=0.90)[1]
Disease.module<-all.genes[1:n.top,]
#
# Subset gene matrix to include only top genes and save it
#
index<-which(rownames(gene.matrix)%in%Disease.module$name)
Disease.module.matrix<-gene.matrix[index,index]
write.table(gene.matrix,file="PF_output/Disease_module_matrix.tsv",row.names=T,
            col.names=T,sep="\t",quote=F)
#
#-------------------------------------------------------------------------------
# Build the graph associated with the top ranking genes and plot it
#-------------------------------------------------------------------------------
#
# Build graph for Disease disease module
#
graph<-graph_from_adjacency_matrix(Disease.module.matrix,mode="undirected",weighted=TRUE)
#
# Color the genes according to the corresponding Expanded Gene List (EGL)
#
Disease.module<-Disease.module[match(rownames(Disease.module.matrix),Disease.module$name), ] # correct the order!
colors<-hcl.colors(nrow(list.matrix),palette="Dark 3",alpha=1)
vertex_colors<-c()
for (i in 1:nrow(Disease.module)) {
  index<-which(rownames(list.matrix)==Disease.module$list.name[i])
  if (length(index)>0) {
    vertex_colors[i]<-colors[index]
  }
}
index<-which(is.na(vertex_colors))
vertex_colors[index]<-"white"
#
# Set the color of nodel labels, to highlight seed genes
#
node.col<-rep("black",nrow(Disease.module))
for (i in 1:nrow(Disease.module)) {
  if("seed"%in%strsplit(Disease.module$source[i],"/")[[1]]) {
    node.col[i]<-"red"
  }
}
#
# Plot the image
#
tiff("PF_output/Disease_module_graph.tiff",width=10,height=10,units="in",res=600,compression="lzw")
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
#
# Save graph for cytoscape
#
edges<-as.data.frame(as_edgelist(graph))
weights<-E(graph)$weight
edges_df<-data.frame(source=edges[,1],target=edges[,2],interaction="interacts_with",
                     weight=weights)
file_name="PF_output/Disease_module_cytoscape.tsv"
write.table(edges_df,file=file_name,sep="\t",row.names=FALSE,quote=FALSE)
#
#-------------------------------------------------------------------------------
# Build a graph with the expanded gene lists, but only for top genes
#-------------------------------------------------------------------------------
#
Disease.module.list.matrix<-ListMatrix(Disease.module)
#
# For each overlap different from zero, calculate a p value by hypergeometric test
#
NB<-nrow(STRING.names) # total background (white and black balls: m+n)
NL<-nrow(Disease.module.list.matrix)
Disease.module.list.matrix.pvalue<-Disease.module.list.matrix
for (i in 1:(NL-1)) {
  for (j in (i+1):NL) {
    Disease.module.list.matrix[i,j]<-as.numeric(Disease.module.list.matrix[i,j])
    if (Disease.module.list.matrix[i,j]==0) {
      Disease.module.list.matrix.pvalue[i,j]<-1
    } else {
      N1<-Disease.module.list.matrix[i,i] # size of EPF_i (number of white balls, m)
      N2<-Disease.module.list.matrix[j,j] # size of EPF_j (number of extracted balls, k)
      NR<-Disease.module.list.matrix[i,j] # size of overlap (number of white balls extracted, q) 
      Disease.module.list.matrix.pvalue[i,j]<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F)) # P(X>NR-1)=P(X>=NR)
      Disease.module.list.matrix.pvalue[i,j]<-as.numeric(formatC(Disease.module.list.matrix.pvalue[i,j],format="e",digits=2))
    }
    Disease.module.list.matrix.pvalue[j,i]<-Disease.module.list.matrix.pvalue[i,j]
  }
}
for (i in 1:NL) {
  Disease.module.list.matrix.pvalue[i,i]<-NA
}
#
# Generate graph from Disease.module.list.matrix
#
temp.matrix<-Disease.module.list.matrix
diag(temp.matrix)<-0 # this is necessary to avoid loops
graph<-graph_from_adjacency_matrix(temp.matrix,mode="undirected",weighted=T)
#
# Generate graph from Disease.module.list.matrix.pvalue
#
temp.matrix<-Disease.module.list.matrix.pvalue
diag(temp.matrix)<-0 # this is necessary to avoid loops
graph.pvalue<-graph_from_adjacency_matrix(temp.matrix,mode="undirected",weighted=T)
#
# Prepare colors for the EGLs
#
colors<-hcl.colors(nrow(Disease.module.list.matrix),palette="Dark 3",alpha=1)
#
# Prepare width for edges
#
edge.w<-E(graph)$weight
edge.w.scaled<-(edge.w-min(edge.w))/(max(edge.w)-min(edge.w))*8+5
#
# Plot it
#
tiff("PF_output/Disease_module_list_matrix.tiff",width=10,height=10,units="in",res=600,compression="lzw")
layout <- layout_with_kk(graph, weights = E(graph)$weight)
plot(graph,
     layout=layout,          
     vertex.size=30,
     vertex.label.cex=2.5,
     vertex.label=paste0(rownames(Disease.module.list.matrix),"\n(size = ",diag(Disease.module.list.matrix),")"),
     vertex.frame.color="black",
     vertex.color=colors,
     edge.label=paste0("overlap: ",E(graph)$weight,"\n p = ",E(graph.pvalue)$weight),
     edge.label.cex=2.0,
     edge.width=edge.w.scaled,
     asp=1,
)
dev.off()
#
#-------------------------------------------------------------------------------
# Gene-set analysis (GSEA) 
#-------------------------------------------------------------------------------
#
top.results<-10 # top results per database
list.number<-1 # EGLs required per pathway
#
score.type<-"score.RWR" # select the score you want to employ for gene-set enrichment analysis
GSEA.RWR<-GSEA.fun(all.genes,top.results,list.number,score.type)
file.name<-"PF_output/GSEA/GSEA_RWR_raw.tsv"
write.table(GSEA.RWR,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# Over-representation analysis (ORA) 
#-------------------------------------------------------------------------------
#
top.results<-10 # top results per database
list.number<-1 # EGLs rquired per pathway
ORA.RWR<-ORA.fun(Disease.module,top.results,list.number)
#
file.name<-"PF_output/ORA/ORA_RWR.tsv"
write.table(ORA.RWR,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
#
#-------------------------------------------------------------------------------
# ORA visualization
#-------------------------------------------------------------------------------
#
for (ORA.type in c("hsa","GO","R-HSA","DO")) {
  myORA<-subset.data.frame(ORA.RWR,grepl(ORA.type,ORA.RWR$ID))
  terms<-myORA[,c("ID","Description","gene.name","p.adjust")]
  terms$Category<-rep("BP",nrow(terms))
  colnames(terms)<-c("ID","Term","Genes","adj_pval","Category")
  terms<-subset.data.frame(terms,select=c("Category","ID","Term","Genes","adj_pval"))
  #
  for (i in 1:nrow(terms)) {
    terms$Genes[i]<-paste0(str_split(terms$Genes[i],"/")[[1]],collapse=",")
  }
  all_names<-unique(unlist(strsplit(paste0(terms$Genes,collapse=","),",")))
  gene_expression<-all.genes[which(all.genes$name%in%all_names),]
  gene_expression<-subset.data.frame(gene_expression,select=c("name","score.RWR"))
  colnames(gene_expression)<-c("ID","logFC")
  gene_expression$logFC<-as.numeric(gene_expression$logFC)
  circ<-circle_dat(terms,gene_expression)
  chord<-chord_dat(circ,gene_expression,terms$Term[1:3])
  if (ORA.type=="hsa") {
    title.str<-"PF_output/ORA/KEGG enrichment"
  } else if (ORA.type=="GO") {
    title.str<-"PF_output/ORA/GO enrichment"
  } else if (ORA.type=="R-HSA") {
    title.str<-"PF_output/ORA/Reactome enrichment"
  } else {
    title.str<-"PF_output/ORA/DO enrichment"
  }
  tiff(paste0(title.str,".tiff"),width=35,height=35,units="cm",res=600,compression="lzw")
  p<-GOChord(chord,space=0.001,title=title.str, 
             gene.order="logFC",
             gene.size = 5,
             process.label = 10,
             lfc.col=c('red','orange','yellow'))
  print(p)
  dev.off()
}
#
#-------------------------------------------------------------------------------
# ORA Pathway visualization (KEGG)
#-------------------------------------------------------------------------------
#
pathway<-subset.data.frame(ORA.RWR,grepl("hsa",ORA.RWR$ID))$ID[1:10]
#
for (i in 1:length(pathway)) {
  p<-pathview(gene.data=as.character(Disease.module$NCBI.id),
              pathway.id=pathway[i],
              species="hsa",
              out.suffix="Disease_ORA",
              gene.idtype="entrez",
              kegg.native=T,
              node.sum="max",
              bins=list(gene=2,cpd=2),
              res=300,
              width=3000,
              height=2500)
}
#
# Move to desired directory
#
files<-list.files(pattern="^hsa.*",full.names=T)
for (i in 1:length(files)) {
  if (grepl("Disease_ORA",files[i])) {
    file.rename(files[i],file.path("PF_output/ORA/KEGG",files[i]))   
  } else {
    file.remove(files[i])
  }
}
#
#-------------------------------------------------------------------------------
# ORA Pathway visualization (Reactome)
#-------------------------------------------------------------------------------
#
pathway<-subset.data.frame(ORA.RWR,grepl("R-HSA",ORA.RWR$ID))$ID
for (i in 1:length(pathway)) {
  entrez_ids<-subset.data.frame(ORA.RWR,ID==pathway[i])$geneID
  entrez_ids<-str_split(entrez_ids,"/")[[1]]
  writeLines(entrez_ids,paste0("PF_output/ORA/Reactome/",pathway[i],".txt"))
}
#
# submit this file at https://reactome.org/PathwayBrowser/#TOOL=AT
#
#-------------------------------------------------------------------------------
# ORA for tissue-specific gene expression
#-------------------------------------------------------------------------------
#
Tissue.ORA(Disease.module)

