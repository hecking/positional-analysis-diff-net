require(igraph)
require(magrittr)
#========================================================
# Remove bots
#------------
# Posts that reacted to a post faster than minLatency
# and that are not referenced themselves
#========================================================
removeBots <- function(g, minLatency) {
  V(g)$indegree <- degree(g,V(g),mode="in")
  V(g)$outdegree <- degree(g,V(g),mode="out")
  sinks <- V(g)[indegree == 1 & outdegree == 0]
  potential_bots <- ends(g,E(g)[latencyMinutes < minLatency])[,2] # posts with an extremely fast response time
  toBeDeleted <- intersect(potential_bots, sinks$name)
  g <- delete_vertices(g, V(g)[is.element(name, toBeDeleted)])
  g
}

#========================================================
# Merge Revisions of Homepages to one instance
#========================================================
mergeRevisions <- function(g, includeTimestamp=FALSE) {
  print("merge revisions")
  if (includeTimestamp) {
    mapping <- as.numeric(as.factor(paste(V(g)$title, V(g)$timemillis)))  
  } else {    
    mapping <- as.numeric(as.factor(V(g)$title))
  }
  g <- contract.vertices(g,mapping,vertex.attr.comb=list(timemillis="min", "first"))
  
  if(!is_dag(g)){
    endpoints <- ends(g, E(g))
    toDelete <- as.numeric(V(g)[endpoints[,1]]$timemillis) > as.numeric(V(g)[endpoints[,2]]$timemillis)
    g <- delete_edges(g, which(toDelete))
  }
  simplify(g)
}

#========================================================
# repair timestamps
#========================================================
repairTimestamps <- function(g) {
  endpoints <- ends(g, E(g))
  negativeEdges <- (as.numeric(V(g)[endpoints[,2]]$timemillis) - 
                      as.numeric(V(g)[endpoints[,1]]$timemillis)) < 0
  
  problematicSources <- unique(endpoints[negativeEdges, 1])
  newTS <- sapply(problematicSources, function(src) {
    neighbourTS <- V(g)[endpoints[endpoints[,1] == src, 2]]$timemillis
    min(neighbourTS) 
  })
  
  V(g)[problematicSources]$timemillis <- newTS
  
  g
}

#========================================================
# merge cycles
#========================================================
mergeCycles <- function(g) {
  #while(!is.dag(g)) {
  
  cycles <- graph.get.subisomorphisms.vf2(g, graph.ring(2, directed=TRUE))
  while(length(cycles) > 0) {
    V(g)$name[cycles[[1]]] <- cycles[[1]][1]
    mapping <- as.numeric(as.factor(V(g)$name))
    g <- contract.vertices(g, mapping, vertex.attr.comb = "first")
    g <- simplify(g)
    cycles <- graph.get.subisomorphisms.vf2(g, graph.ring(2, directed=TRUE))
  }
  
  cycles <- graph.get.subisomorphisms.vf2(g, graph.ring(3, directed=TRUE))
  while(length(cycles) > 0) {
    
    V(g)$name[cycles[[1]]] <- cycles[[1]][1]
    mapping <- as.numeric(as.factor(V(g)$name))
    
    g <- contract(g, mapping, vertex.attr.comb = "first")
    g <- simplify(g)
    cycles <- graph.get.subisomorphisms.vf2(g, graph.ring(3, directed=TRUE))
  }
  
  cycles <- graph.get.subisomorphisms.vf2(g, graph.ring(4, directed=TRUE))
  while(length(cycles) > 0) {
    
    V(g)$name[cycles[[1]]] <- cycles[[1]][1]
    mapping <- as.numeric(as.factor(V(g)$name))
    
    g <- contract(g, mapping, vertex.attr.comb = "first")
    g <- simplify(g)
    cycles <- graph.get.subisomorphisms.vf2(g, graph.ring(4, directed=TRUE))
  }
  
  #  g <- simplify(g)
  #}
  
  g
}

#========================================================
# count and print how many nodes exist of any media type
#========================================================
printNumberNodeTypes <- function(g,minLatency,merge,botRemoval) {
  tweets <- V(g)[type=="TWEET"]
  retweets <- V(g)[type=="RETWEET"]
  webpages <- V(g)[type=="WEB"]
  wikiArticles <- V(g)[type=="WIKIPEDIA"]
  youtubePosts <- V(g)[type=="YOUTUBE"]
  facebookPosts <- V(g)[type=="FACEBOOK"]
  
  print(paste("Tweets", length(tweets)))            
  print(paste("Retweets", length(retweets)))            
  print(paste("Webpages", length(webpages)))            
  print(paste("Wikipedia", length(wikiArticles)))            
  print(paste("Youtube", length(youtubePosts)))            
  print(paste("Facebook", length(facebookPosts)))            
}


#========================================================
# cleanData
#========================================================
cleanData <- function(g,minLatency,merge,botRemoval) {
  print("Nodes before cleaning")
  printNumberNodeTypes(g)
  
  g <- g %>% repairTimestamps %>% mergeRevisions(TRUE) %>% mergeCycles %>% weightEdgesTemporally
  if(merge){
    g <- g %>% mergeRevisions %>% weightEdgesTemporally
  }
  if(botRemoval){
    g <- removeBots(g,minLatency)
  }
  # delete isolated nodes
  V(g)$degree <- degree(g,V(g))
  g <- delete_vertices(g, V(g)[degree == 0])
  if (!is.dag(g)) {
    stop("The graph is not a DAG. Further cleaning necessary")
  }
  
  print("-----")
  print("Nodes after cleaning")
  printNumberNodeTypes(g)
  
  g
}