import os as os
def  nameshift(pre,ext,si,ei,offset):
  """  
  shift the file name number
  Input parameter :
  pre    : prefix of the filename
  ext    : extension of the filename
  si     : starting index
  ei     : end index
  offset : offset number
  """


  if si>ei:
    ss=si
    si=ei
    ei=ss

  if offset==0: 
    return
  elif offset>0:
    if (ei+offset)>999: 
      return
    else:
      ii=range(ei,si-1,-1)    # large number first
  else:
    if (si+offset)<0: 
        return
    else:
        ii=range(si,ei+1)     # small number first
    
  for i in ii:
    fname=pre+'%03d'%i+'.'+ext
    fname1=pre+'%03d'%(i+offset)+'.'+ext
    os.rename(fname,fname1)
