library('stringr')

mask = c(0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,
         1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1)

base = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

enabled = length(which(mask==1))

outdir = '/mnt/raid6_data/user/aranyics/rdcm2/vdmodel_4node_s01/'

for ( i in 0:(2^enabled-1) )
{
  bitvect = base
  bitvect[which(mask==1)] = rev(as.numeric(intToBits(i))[1:enabled])
  #if (i == 20) break
  write.table(str_c(bitvect[1:16], collapse=''), str_c(outdir, 'A.csv'), append=TRUE, row.names=FALSE, col.names=FALSE)
  write.table(str_c(bitvect[17:64], collapse=''), str_c(outdir, 'B.csv'), append=TRUE, row.names=FALSE, col.names=FALSE)
}

