library(GO.db)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
head(keys(org.Hs.eg.db))
keytypes(GO.db)
keytypes(org.Hs.eg.db)

   #------------------------------------------------------------
   # get all entrezid gene go annotations, all 3 ontologies
   #------------------------------------------------------------

all.entrez <- keys(org.Hs.eg.db)
tbl.gene <- select(org.Hs.eg.db, keys=all.entrez, columns=c("SYMBOL", "GO", "ENTREZID"), keytype="ENTREZID")
dim(tbl.gene)  # 379,756 5

   #------------------------------------------------------------
   # add ensembl ids
   #------------------------------------------------------------

length(tbl.gene$ENTREZID) # 378756
length(unique(tbl.gene$ENTREZID)) # 61547


tbl.entrezEnsembl <- select(EnsDb.Hsapiens.v79, key=(unique(tbl.gene$ENTREZID)),
                            keytype="ENTREZID",
                            columns=c("GENEID"))
dim(tbl.entrezEnsembl)  # 28257 2
dups <- which(duplicated(tbl.entrezEnsembl$ENTREZID))
length(dups)  # 4024   are any of these worthy?  let's assume not for now
tbl.entrezEnsembl <- tbl.entrezEnsembl[-dups,]
dim(tbl.entrezEnsembl)  # 24233 2

   #----------------------------------------------------------------
   # now merge tbl.gene, with entrez and go, to include ensg column
   #----------------------------------------------------------------

tbl.geneGo <- merge(tbl.gene, tbl.entrezEnsembl, by="ENTREZID")
dim(tbl.geneGo)  # 335712 6
dim(unique(tbl.geneGo))  # 335712 6

tbl.geneGoOut <- tbl.geneGo[, c("GENEID", "GO", "ONTOLOGY", "EVIDENCE")]
tbl.geneGoOut$interaction <- sprintf("GO_%s", tbl.geneGo$ONTOLOGY)

write.table(tbl.geneGoOut,
            sep="\t",
            file="../neo4j/import/ensgGO.tsv",
            quote=FALSE,
            row.names=FALSE)

   #----------------------------------------------------------------
   # associate goid with goTerm, write these out as nodes
   #----------------------------------------------------------------

tbl.goTerms <- unique(select(GO.db, keys=tbl.geneGo$GO, keytype="GOID", columns=c("TERM", "ONTOLOGY")))

dim(tbl.geneGo)    # 335712 6
dim(tbl.goTerms)   # 18309 3
tbl.goTerms$type = "GO"

write.table(tbl.goTerms,
            sep="\t",
            file="../neo4j/import/goidTermsOntology.tsv",
            quote=FALSE,
            row.names=FALSE)



   #----------------------------------------------------------------
   # add the goterm column to the gene table
   #----------------------------------------------------------------

tbl.geneGoWithTerms <- merge(tbl.geneGo, tbl.goTerms[, 1:2], by.x="GO", by.y="GOID")
dim(tbl.geneGoWithTerms)  # 335712 7
table(tbl.geneGoWithTerms$ONTOLOGY)
#     BP     CC     MF
# 155759  96196  79015

   #----------------------------------------------------------------
   #  extract and write out the GO terms, ready for import
   #----------------------------------------------------------------

tbl.nodes.go <- unique(tbl.geneGoWithTerms[, c("GO", "ONTOLOGY", "TERM")])
dim(tbl.nodes.go)  # 18309 3
table(tbl.nodes.go$ONTOLOGY)
   #    BP    CC    MF
   # 12315  1754  4239
write.table(tbl.nodes.go, file="tbl.nodes.go.tsv", sep="\t", quote=FALSE, row.names=FALSE)

   #----------------------------------------------------------------
   #  need a table of ensg, symbol, uniprot
   #----------------------------------------------------------------

tbl.nodes.genes <- select(EnsDb.Hsapiens.v79, key=keys(EnsDb.Hsapiens.v79),
                          keytype="GENEID", columns=c("GENEID", "SYMBOL", "UNIPROTID"))
deleters <- which(is.na(tbl.nodes.genes$UNIPROTID))
length(deleters)
tbl.nodes.genes <- unique(tbl.nodes.genes[-deleters,])
dim(tbl.nodes.genes)  # 74485 3

   #----------------------------------------------------------------
   # write out three sets of nodes, and two sets of edges:
   #   tbl.nodes.ensg
   #   tbl.nodes.symbol
   #   tbl.nodes.uniprot
   #   tbl.edges.ensg.symbol
   #   tbl.edges.ensg.uniprot
   #----------------------------------------------------------------


tbl.ensgSymbols <- unique(tbl.nodes.genes[, c("GENEID", "SYMBOL")])
tbl.ensgSymbols$type <- "ensg"
write.table(tbl.ensgSymbols, file="../neo4j/import/ensgSymbols.tsv", quote=FALSE, row.names=FALSE, sep="\t")

tbl.uniprotSymbols <- unique(tbl.nodes.genes[, c("UNIPROTID", "SYMBOL")])
tbl.uniprotSymbols$type <- "uniprot"
dim(tbl.uniprotSymbols)
write.table(tbl.ensgSymbols, file="../neo4j/import/uniprotSymbols.tsv", quote=FALSE, row.names=FALSE, sep="\t")

xxxx





tbl.nodes.ensg <- unique(tbl.nodes.genes[, 1:2])
dim(tbl.nodes.ensg)  # 65774 2
tbl.nodes.uniprot <- unique(tbl.nodes.genes[, 3, drop=FALSE])
dim(tbl.nodes.uniprot) # 70196 1

tbl.edges.ensg.uniprot <- unique(tbl.nodes.genes[, c("GENEID", "UNIPROTID")])
dim(tbl.edges.ensg.uniprot)  # 132907
deleters <- which(is.na(tbl.edges.ensg.uniprot$UNIPROTID))
length(deleters)
tbl.edges.ensg.uniprot <- tbl.edges.ensg.uniprot[-deleters,]
dim(tbl.edges.ensg.uniprot)  # 74484 2

tbl.edges.ensg.go <- unique(tbl.geneGoWithTerms.


length(all.entrez)  # 61547


dim(tbl.entrezEnsembl) # 20661 2


tbl.go <- select(GO.db, keys=unique(tbl.gene$GO), keytype="GOID", columns=c("TERM", "DEFINITION"))
dim(tbl.go)
tbl.merged <- unique(merge(tbl.gene, tbl.go, by.x="GO", by.y="GOID"))
dim(tbl.merged) # [1] 180,038 8
write.table(tbl.merged, file="tbl.geneGO.tsv", sep="\t", quote=FALSE, row.names=FALSE)
system("wc -l tbl.geneGO.tsv") # 180k

tbl.nodes <- data.frame(id=unique(c(tbl.merged$GO, tbl.merged$ENSEMBL)), type="unassigned", stringsAsFactors=FALSE)
dim(tbl.nodes)  # 33424
type <- tbl.nodes$type
type[grep("GO:", tbl.nodes$id)] <- "GOBP"
type[grep("ENSG", tbl.nodes$id)] <- "gene"
tbl.nodes$type <- type
head(tbl.nodes)
table(tbl.nodes$type)
subset(tbl.nodes, type=="unassigned")
dim(tbl.nodes)
tbl.nodes[12490:12510,]
tbl.nodes <- subset(tbl.nodes, type != "unassigned")

head(tbl.merged)

tbl.edges <- tbl.merged [, c("ENSEMBL", "GO", "EVIDENCE", "SYMBOL")]
dim(tbl.edges)  # 180038
all(tbl.edges$GO %in% tbl.nodes$id)
all(tbl.edges$ENSEMBL %in% tbl.nodes$id)
setdiff(tbl.edges$ENSEMBL, tbl.nodes$id)
dim(subset(tbl.edges, is.na(ENSEMBL)))
head(subset(tbl.merged, is.na(ENSEMBL)))


tbl.ens <- select(EnsDb.Hsapiens.v79, key=(tbl.merged$ENTREZID),
                                keytype="ENTREZID",
                                columns=c("GENEID", "SYMBOL"))
dim(tbl.ens) # 20661
length(grep("ENSG", tbl.ens$GENEID))   # 20106
tbl.ens <- tbl.ens[grep("ENSG", tbl.ens$GENEID),]
dim(tbl.ens)  # 20106 3






