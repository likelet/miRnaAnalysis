library(VennDiagram)
library(RColorBrewer)
input_dir = "."
output = "All_pair_venn.pdf"
file_lst = list.files(input_dir, pattern = "*.tsv$")
n_file = length(file_lst)
target_lst = list()
gene_names = array()
for (i in 1:n_file) {
  mat <- read.table(file_lst[i], header = T)
  target_lst[[i]] = paste(mat[1:2,1], mat[1:2,2], sep = "--")
  names(target_lst)[i] = sub(".tsv", "", file_lst[i])
  gene_names = c(gene_names, levels(mat[,1]))
}
gene_names = unique(na.omit(gene_names))
fig <- venn.diagram(target_lst, fill = brewer.pal(n_file,"Dark2"),
             cex=2, filename=NULL)
pdf(output)
grid.draw(fig)
dev.off()

venn_gene_plot = function(gene){
  miRNA = list()
  for (i in 1:n_file) {
    mat <- read.table(file_lst[i], header = T)
    miRNA[[i]] = mat[which(mat[,1] == gene),2]
    names(miRNA)[i] = sub(".tsv", "", file_lst[i])
  }
  fig <- venn.diagram(miRNA, fill = brewer.pal(n_file,"Dark2"),
                      cex=2, filename=NULL)
  pdf(paste(gene, "_venn.pdf", sep = ""))
  grid.draw(fig)
  dev.off()
}

for (i in gene_names){
  venn_gene_plot(i)
}
