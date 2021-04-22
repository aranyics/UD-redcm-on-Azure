redcm_mask_matrix2id = function(M, mask, digits)
{
  
  if (class(M) == "character")
    M = matrix(as.numeric(strsplit(M, split = '')[[1]]), nrow = sqrt(nchar(M)))
  
  M = as.logical(M) * 1
  
  id = paste( strtoi( paste(M[mask==1], collapse = ''), 2 ) )
  if (id == 'NA') id = '0' #if mask is empty, there are no bits to change
  
  for (i in nchar(id):(digits-1))
    id = paste('0', id, sep = '')
  
  return (id)
  
}

redcm_mask_id2matrix = function(id, mask, Mbase=NULL)
{
  
  M = mask
  n = sqrt(length(mask))
  
  bits = tail(as.numeric(rev(intToBits(as.integer(id)))), length(which(mask==1)))
  
  M[mask == 1] = bits
  M = matrix(M, nrow = n)
  
  if (!is.null(Mbase))
    M = M + Mbase #assuming mask and Mbase are disjunct bitwise
  
  return (M)
  
}

redcm_mask_mid2matrix = function(id.list, mask.list, base.list, digits)
{
 
  if (length(id.list) == 1 && class(id.list) == "character")
  {
    mid = as.character(id.list)
    for( i in nchar(mid):(length(mask.list)*digits-1) )
      mid = paste0('0', mid)
    id.list = as.list(strsplit(gsub(sprintf("(.{%d})", digits), "\\1 ", mid), ' ')[[1]])
  }
  else if (length(id.list) != length(mask.list) || length(id.list) != length(base.list) )
  {
    #wrong input length
    return (NULL)
  }
  
  Ml = list()
  
  for (i in 1:length(id.list))
  {
    Ml[[i]] = redcm_mask_id2matrix( id.list[[i]], mask.list[[i]], base.list[[i]] )
  }
  
  return(Ml)
  
}
