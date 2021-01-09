library(GO.db)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)

tbl.geneGo <- select(org.Hs.eg.db, keys="DHODH", columns=c("SYMBOL", "GO", "ENTREZID"), keytype="SYMBOL")
tbl.goTerms <- unique(select(GO.db, keys=tbl.geneGo$GO, keytype="GOID", columns=c("TERM", "DEFINITION")))
