#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (subject)", call.=FALSE)
}

## -----------------------------------------------------------------------------------------
## set environmental variables
.libPaths(rev(dir('~/R/x86_64-pc-linux-gnu-library', full.names=TRUE)))
cat( 'R package path: ' )
cat( .libPaths(), '\n', sep=':')

suppressPackageStartupMessages( library("ReDCM") )

subject = args[1]
Amask = c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0)
Abase = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
Bmask = t( array(
          c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1),
        dim = c(16, 3)) )
digits = 3

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

denan.mx = function(mx)
{
  mx[is.nan(mx)] = 0
  return (mx)
}

source(paste(dirname(thisFile()), '/redcm_mask_mid.R', sep = ''))

dcmdir = paste(dirname( thisFile() ), '/../../../results/modelspace/vdmodel_4node_', subject, '/dcm/', sep = '')
outdir = paste(dcmdir, '../dcm_ms/', sep='')
dir.create(outdir, showWarnings = FALSE)
d = list.files(dcmdir, recursive = TRUE)

for (i in 1:length(d))
{
  id = tail( strsplit(d[i], '_')[[1]], 1 )
  #cat(d[i])
  load( paste(dcmdir, d[i], sep = '') )
  
  if (i == 1)
  {
    nn = DCMe@M[[1]]@l
    nm = DCMe@M[[1]]@m
    cols = c(apply(array(1:nn), 1, function(x){rep(x, nn)}))
    rows = rep(1:nn,nn)
    Bs = apply(array(1:nm), 1, function(x){paste('B', x, sep = '')})
    Cs = apply(array(1:nm), 1, function(x){paste('C', x, sep = '')})
    dcm.hd = c('ID', 'mID', 'A', Bs, 'C', 'Fe')
    dcm.hd = c(dcm.hd, paste(rep('A_', nn*nn), cols, rows, sep=''), paste(rep('pA_', nn*nn), cols, rows, sep=''))
    for (k in 1:nm)
    {
      dcm.hd = c(dcm.hd, paste(rep(paste(Bs[k],'_',sep =''), nn*nn), cols, rows, sep=''), paste(rep(paste('p',Bs[k],'_',sep =''), nn*nn), cols, rows, sep=''))
    }
    for (k in 1:nm)
    {
      dcm.hd = c(dcm.hd, paste(rep(paste(Cs[k],'_',sep =''), nn), 1:nn, sep=''), paste(rep(paste('p',Cs[k],'_',sep =''), nn), 1:nn, sep=''))
    }
    dcm.length = length(dcm.hd)
    dcm.hd = paste(dcm.hd, collapse = ',')
    fileConn<-file(paste(outdir,"fullspace_dcm.csv", sep = ''), "wt")
    writeLines(dcm.hd, fileConn)
    close(fileConn)
  }
  
  aid = redcm_mask_matrix2id(DCMe@Ep@A, Amask, digits)
  bid = apply(array(1:DCMe@M[[1]]@m), 1, function(k) {redcm_mask_matrix2id(DCMe@Ep@B[,,k], Bmask[k,], digits)})
  cid = "000"
  
  dcm = c(id, paste(c(aid, bid), collapse = ''), aid, bid, cid)
  dcm = c(dcm, round(DCMe@Fe, 4))
  dcm = c( dcm, c(round(DCMe@Ep@A, 4)), c(denan.mx(round(DCMe@Pp@A, 4))) )
  for (m in 1:DCMe@M[[1]]@m)
    dcm = c( dcm, c(round(DCMe@Ep@B[,,m], 4)), c(denan.mx(round(DCMe@Pp@B[,,m], 4))) )
  dcm = c( dcm, c(round(DCMe@Ep@C, 4)), c(denan.mx(round(DCMe@Pp@C, 4))) )
  dcm = paste(dcm, collapse = ',')
  fileConn<-file(paste(outdir,"fullspace_dcm.csv", sep = ''), "at")
  writeLines(dcm, fileConn)
  close(fileConn)
  
  if (i %% 1000 == 0)
    cat('-')
  
}

df = read.table(paste(outdir,"fullspace_dcm.csv", sep = ''), sep = ',', header = TRUE, colClasses = 'character')
df = df[order(as.numeric(df$ID)),]
df[,-c(2:7)] = lapply(df[,-c(2:7)], as.numeric)
write.table(df, paste(outdir,"fullspace_dcm.csv", sep = ''), sep = ',', col.names = TRUE, row.names = FALSE)

cat('\n')
