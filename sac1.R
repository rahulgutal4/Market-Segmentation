#install.packages("igraph")
#install.packages("hash")

library(igraph)
library(LSAfun)

cosineSimilarity = function(attr, v){
  cos = matrix(nrow = v, ncol = v)
  for(i in 1:v){
    vec = as.numeric(attr[i,])
    for(j in i:v){
      cos[i,j] = cosine(as.numeric(attr[j, ]), vec)  
      cos[j,i] = cos[i,j]
    }
  }
  return(cos)
}

phase1 = function(gg, alpha, cos){
  vertices = V(gg)
  v = length(vertices)
  c = c(1:v)

  attr = read.table("data/fb_caltech_small_attrlist.csv", header = TRUE, sep = ',')
  # phase 1
  for(p in 1:15){
    
    flag = FALSE
    for(i in 1:v){
      
      sMod1 = modularity(gg, c)
      finalC = 0
      maxDelta = 0
      com = setdiff(unique(c), c[i])
      #com = unique(c)
      for(j in com){
        #cat(p,' ',' ',i,' ',j, '\n')
        
        temp = c[i]
        c[i] = j
        sMod2 = modularity(gg, c)
        aMod = 0.0
        indices = setdiff(which(c == j),i)
        for(ind in indices){
          aMod = aMod + cos[i,ind]  
        }
        
        delta = alpha*(sMod2 - sMod1) + (1-alpha)*aMod
        if(delta > maxDelta){
          maxDelta = delta
          finalC = j
        }
        c[i] = temp
        
      }
      
      if(maxDelta != 0){
        c[i] = finalC
        flag = TRUE
      }
    }
    if(!flag)
      break
  }
  return(c)
}

sac1 = function(alpha = 0){
  
  gg = read_graph("data/fb_caltech_small_edgelist.txt",format= c("edgelist"))
  attr = read.csv("data/fb_caltech_small_attrlist.csv", header = TRUE, sep = ',')
  vertices = V(gg)
  v = length(vertices)

  cos = cosineSimilarity(attr, v)
  c = phase1(gg, alpha, cos)
  ggc = contract.vertices(gg, c)
  ggc = simplify(ggc, remove.multiple = TRUE, remove.loops = TRUE)
  com = unique(c)
  result = list()
  i = 1
  for(comm in com){
    cm = which(c == comm)
    temp = c()
    for(ind in cm){
      if(length(temp) == 0) temp = c(as.numeric(attr[ind,]))
      else temp = temp + as.numeric(attr[ind,])
    }
    result[[i]] = temp
    i = i + 1
  }
  
  v = length(com)
  cos = matrix(nrow = v, ncol = v)
  for(i in 1:v){
    for(j in i:v){
      cos[i,j] = cosine(result[[j]], result[[i]])  
      cos[j,i] = cos[i,j]
    }
  }

  f = file("communities.txt","w")
  com = unique(c)
  for(comm in com){
    cm = which(c == comm)
    cat(as.character(cm), file = f)
    cat("\n", file=f)
  }
  close(f)
}