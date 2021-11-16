# This is moc_gapbk algorithm copied and pasted from the now defunct CRAN
# repostory site (https://cran.r-project.org/web/packages/moc.gapbk/index.html).
# It was used to test why it was taking so long when using its LS methods. The
# results were that it was calculating and evaluating a lot of new solutions
# with its LS methods, since its not bounded by max evaluations or time.

#Para evitar que de una NOTE en pathrelinking con variable i
utils::globalVariables(c("i"))


generate.initial.population<-function(num_objects, num_k, pop_size){
  population.P = t(sapply(1:pop_size, function(x)  sample(1:num_objects, num_k ,replace=F)  ))
  rm(.Random.seed, envir=globalenv())
  return(as.matrix(population.P))
}



generate.groups<-function(pop_size, medoids, dmatrix1, dmatrix2){
  
  groups<-list(1:pop_size)
  formacion<-rep(NA,nrow(dmatrix1))
  
  for (p in 1:pop_size) {
    for (i in 1:nrow(dmatrix1)) {
      grupo<-0
      gen<-i
      valor<-(dmatrix1[gen, medoids[p,]])
      grupo<-unname(which.min(valor))
      formacion[i]<-grupo
    }
    groups[[p]] <- formacion
    names(groups[[p]])<-rownames(dmatrix1)
  }
  
  return(groups)
  
}




verify.feasibility<-function(population, par_destino, num_objects){
  
  #recorre cada individuo, si encuentra un cromosoma repetido(medoide) lo reemplaza aleatoriamente por otro.
  for (p in 1:nrow(population)) {
    
    while(length(which(duplicated(population[p,])))>0){
      id_rep<-which(duplicated(population[p,]))
      num_rep<-length(id_rep)
      
      for (i in 1:num_rep) {
        population[p,id_rep[i]]=sample(1:num_objects,1, replace=F)
      }
      
    }
    
  }
  
  par_destino<-population
  
  return(par_destino)
  
}

singletons.delete<-function(groups_formados, poblacion, num_k) {
  
  i=1
  solucion_singleton<-0
  
  for (individuo in 1:length(groups_formados)) {
    
    for (cluster in 1:num_k) {
      
      if(length(which(groups_formados[[individuo]] == cluster))<2){
        solucion_singleton[i]<-individuo
        i<-i+1
        break
      }
      
    }
    
  }
  
  if(length(solucion_singleton)>0 && min(solucion_singleton) >0){
    groups_formados<-groups_formados[-c(solucion_singleton)]
    poblacion<-poblacion[-c(solucion_singleton),]
    
  }
  
  
  return(list(groups=groups_formados, poblacion=poblacion))
  
  
}





singletons.repair<-function(population.repair, dmatrix1, dmatrix2, num_objects) {
  
  
  table.repair.groups<-generate.groups(nrow(population.repair), population.repair, dmatrix1, dmatrix2)
  
  #si encuentra solucion con cluster singletons, quita el medoide que produce eso y lo reemplazo por uno aleatorio.
  #el reemplazo se verifica para evitar que sea un medoide ya incluido y ademas que no produzca otros singletons
  
  genes_total <- c(1:num_objects)
  
  for (individuo in 1:nrow(population.repair)) {
    
    #puede suceder que alterar por ejemplo el ultimo medoide del cromosoma hago que el medoide 1 del mismo
    #cromosoma ahora produzca singletons. Por eso cada vez que hay singletons vuelve a verificar desde el primer medoide del cromosoma
    columna <- 1
    
    while(columna <= ncol(population.repair)){
      
      todoOK=FALSE
      
      while(todoOK==FALSE){
        
        if(length(which(unlist(table.repair.groups[[individuo]]) == columna))<3){
          
          #cat("individuo ", individuo, " medoide ", columna, " produce singletons", "\n")
          #obtiene los medoides del individuo donde se encontraron singletons
          genes_cromosoma <- population.repair[individuo,1:ncol(population.repair)]
          #identifica cuales de todo el conjunto de genes no estan en el cromosoma,
          #los cuales serviran para el reemplazo del medoide que es singleton
          posibles_reemplazos<-setdiff(genes_total,intersect(genes_total, genes_cromosoma))
          #de los genes reciente identificados escojo uno aleatoriamente
          reemplazo<-sample(posibles_reemplazos,1, replace=F)
          
          #reemplazo tiene el gene que sera ahora usado como medoide
          population.repair[individuo,columna]=reemplazo
          #re forma los groups con el nuevo medoide.
          table.repair.groups[[individuo]]<-generate.groups(1, as.matrix(t(population.repair[individuo,])), dmatrix1, dmatrix2)
          todoOK=FALSE
          columna=1
          
        }else{
          todoOK=TRUE
          columna=columna+1
        }
        
      }
      
      
    }
    
    #el individuo i debe estar mejorado cuando haya salido de todos los while, o sea
    #el population.repair[individuo, ] y tablatable.repair.groups[[individuo]] debe estar con dos genes en cada grupo
  }
  
  #groups_sin_singletons<-subset(groups_formados, !(groups_formados %in% cluster_singleton))
  
  
  
  return(list("integer" = table.repair.groups, "integer"=population.repair))
  
}





generate.crossover.k.points<-function(pop_size, num_k, population.from, population.to, rat_cross){
  
  p=0.50
  
  for (i in seq(1, pop_size, by = 2)) {
    
    parejas<-sample(1:pop_size,2, replace=F)#SI LO CAMBIO a TRUE puede ocasionar parejas iguales, o sea 3, 3 por ejemplo
    
    cruce=stats::runif(1,0,1)
    
    if (cruce < rat_cross) {
      
      for (j in 1:num_k) {
        
        aleatorio=stats::runif(1,0,1)
        if (aleatorio<=p) {
          population.to[i,j]=as.matrix(population.from[parejas[1],j])
          population.to[(i+1),j]=as.matrix(population.from[parejas[2],j])
        }else{
          population.to[i,j]=as.matrix(population.from[parejas[2],j])
          population.to[(i+1),j]=as.matrix(population.from[parejas[1],j])
        }
      }
      
      
    }else{
      
      population.to[i,1:num_k]<-as.matrix(population.from[parejas[1],1:num_k])
      population.to[(i+1),1:num_k]<-as.matrix(population.from[parejas[2],1:num_k])
      
    }
    
    
  }
  
  return(population.to)
  
}




generate.mutation.random.controller<-function(pop_size, num_k, population, rat_muta, num_objects){
  
  for (p in 1:pop_size) {
    
    muta=stats::runif(1,0,1)
    
    if (muta < rat_muta) {
      posicion_mutar = sample(1:num_k,1, replace=F)
      medoide=sample(1:num_objects,1, replace=F)
      population[p,posicion_mutar]=medoide
    }
    
  }
  
  return(population)
  
}




calculate.objective.functions <- function(pop_size, population, groups, par_distancia, num_k, local_search) {
  
  
  #pop_size=pop_size
  #population=population.P
  #groups=table.groups.P
  #dmatrix1=dmatrix1
  #num_k=num_k
  #local_search=FALSE
  
  
  if (local_search==TRUE) {
    population= population[,2:(num_k+1)]
  }
  
  
  
  fobj<-rep(0,length(groups))
  
  
  
  for (x in 1:length(groups)) {
    
    fitness_counter <<- fitness_counter + 1
    message(str_interp("Fitness counter: ${fitness_counter}"))
    
    ############### INICIALIZA MATRICES """""""
    n= nrow(par_distancia) #Numero de elementos
    K=num_k #Numero de centros
    D <- matrix(0, K, n) #Matriz de distancia centros x elementos
    
    # XB varianza
    
    #Calcula Distancia entre cada medoide y los genes
    for (k in 1:K) {
      for (i in 1:n) {
        medoide_cluster=population[x, k]
        D[k, i] = (par_distancia[medoide_cluster, i])  ^ 2
      }
    }
    
    
    
    XB_numerador = sum( D )
    
    
    # XB separacion
    
    #genera cada par de medoide
    pares_medoides= t(utils::combn(population[x,1:num_k],2))
    
    s=vector()
    for (fila in 1:nrow(pares_medoides)) {
      s[fila]= (par_distancia[pares_medoides[fila,1],pares_medoides[fila,2]]) ^ 2
    }
    separacion=min(s)
    
    
    XB_denominador = separacion
    
    
    
    
    ###varianza=0
    
    ###for (k in 1:num_k) {
    #distancia entre cada elemento y el medoide k
    ###elementos_cluster=unname(which(unlist(groups[[x]])==k))
    ###medoide_cluster=population[x,k]
    #suma las distancias para todos los k
    ###varianza=varianza+sum(par_distancia[medoide_cluster,elementos_cluster]^2)
    ###}
    
    #genera cada par de medoide
    ###pares_medoides= t(combn(population[x,1:num_k],2))
    
    ###s=vector()
    ###for (fila in 1:nrow(pares_medoides)) {
    ###s[fila]=par_distancia[pares_medoides[fila,1],pares_medoides[fila,2]]^2
    ###}
    ###separacion=min(s)
    
    ###fobj[[x]]=varianza/(num_k*separacion)
    
    fobj[[x]]=XB_numerador/(n*XB_denominador)
    
    
  }
  
  
  return(unlist(fobj))
  
}


calculate.ranking.crowding<-function(pop_size, population, groups, number_objectives, local_search, dmatrix1, dmatrix2, num_k){
  
  
  obj1<-calculate.objective.functions(pop_size, population, groups, dmatrix1, num_k, local_search)
  obj2<-calculate.objective.functions(pop_size, population, groups, dmatrix2, num_k, local_search)
  
  population <- cbind(population,obj1,obj2)
  
  if(local_search==TRUE){
    
    if (nrow(population)<=1) {
      objetivos=as.matrix(t(population[,(num_k+2):(num_k+number_objectives+1)]))
    }else{
      objetivos=population[,(num_k+2):(num_k+number_objectives+1)]
    }
    
    
  }else{
    
    if (nrow(population)<=1) {
      objetivos=as.matrix(t(population[,(num_k+1):(num_k+number_objectives)]))
    }else{
      objetivos=population[,(num_k+1):(num_k+number_objectives)]
    }
    
  }
  
  ranking <- nsga2R::fastNonDominatedSorting(objetivos)
  
  paretoranking <- calculate.ranking(pop_size,ranking)
  population <- cbind(population, paretoranking)
  
  #calcula crowding Poblacional en el Frente
  objRange <- calculate.objectives.range(population,num_k,number_objectives,local_search)
  cd <- nsga2R::crowdingDist4frnt(population,ranking,objRange)
  population <- cbind(population,crowding=c(apply(cd,1,sum)))
  population<-as.data.frame(population)
  population<-population[order(population$paretoranking, -(population$crowding)), ]
  population<-as.matrix(population)
  
  #devuelve la poblacion ordenada segun Jerarquia y crowding.
  
  return(population)
  
  #population es un objeto matrix
  
}

calculate.ranking<-function(pop_size, ranking){
  
  rankIndex <-integer(pop_size)
  i <- 1
  while (i <= length(ranking)) {
    rankIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  
  return(rankIndex)
}

calculate.objectives.range<-function(population, par_var, par_obj, local_search){
  
  
  
  if(local_search==TRUE){
    
    if (nrow(population)<=1) {
      objetivos=as.matrix(t(population[,(par_var+2):(par_var+par_obj+1)]))
    }else{
      objetivos=population[,(par_var+2):(par_var+par_obj+1)]
    }
    
    rangoObj<-apply(objetivos, 2, max) -  apply(objetivos, 2, min)
    
  }else{
    
    if (nrow(population)<=1) {
      objetivos=as.matrix(t(population[,(par_var+1):(par_var+par_obj)]))
    }else{
      objetivos=population[,(par_var+1):(par_var+par_obj)]
    }
    
    rangoObj<-apply(objetivos, 2, max) -  apply(objetivos, 2, min)
  }
  
  return(rangoObj)
  
}



generate.pareto.local.search<-function(population_pareto, neighborhood, num_k, num_objects, pop_size, dmatrix1, dmatrix2, number.objectives){
  
  
  #Dejar solo medoides. Pero como maximo el tama?o de la poblaci?n.
  if(nrow(population_pareto)>pop_size){
    poblacion_A0<-population_pareto[1:pop_size,1:num_k]
  }else{
    poblacion_A0<-population_pareto[,1:num_k]
  }
  
  
  #Dejar solo medoides sin repeticion
  poblacion_A0<-poblacion_A0[!duplicated(poblacion_A0[,1:num_k]),]
  
  #Marcar como no exploradas todas las soluciones Pareto. 0=FALSE, 1=TRUE
  poblacion_A0<-cbind(explored=rep(0,nrow(poblacion_A0)),poblacion_A0)
  
  names(poblacion_A0)<-paste("V",1:(num_k+1),sep = "")
  
  #Crear poblacion Archivo A
  poblacion_A<-poblacion_A0
  
  #cat("\n poblacion antes \n")
  #print(poblacion_A)
  
  #Establece tamano vencindad
  num_neighborhood<-ceiling(neighborhood*num_objects)
  #poblacion vecindad
  
  
  #solo se se realiza PLS si la poblacion a mejorar tiene mas de 2 soluciones
  if(nrow(poblacion_A0)>1){
    
    
    
    #poblacion A0 tiene todas soluciones de Pareto, pero se escoge solo 1
    #y sobre esa se hace el procedimiento de PLS
    #Esto se realiza en poblacion_A0<-subset(poblacion_A, V1=='0')
    
    while(nrow(poblacion_A0)>1){
      
      
      posicion_alterar<-sample(1:num_k,1)
      s=sample(1:nrow(poblacion_A0), 1)
      
      solucion_S<-as.data.frame(poblacion_A0[s,])
      
      #almacena los medoides que no estan incluidos en archivo
      #la idea es que al alterar la poblacion A0, no se repitan los medoides
      lista_medoides<-as.vector(as.matrix(poblacion_A0[,2:(num_k+1)]))
      lista_medoides_unicos<-unique(lista_medoides)
      lista_medoides_total<-1:num_objects
      lista_medoides_disponibles<-setdiff(lista_medoides_total,lista_medoides_unicos)
      
      #puedo llegar a tener todos medoides explorados y quedarme sin medoide para asignar como nuevo medoide
      #por eso solo aplico cuando haya al menos un medoide disponible en lista_medoides_disponibles
      ###if(length(lista_medoides_disponibles)>0){
      
      #Lleno la poblacion de neighborhood N
      for (i in 1:num_neighborhood) {
        
        nuevo_medoide<-lista_medoides_disponibles[sample(1:length(lista_medoides_disponibles),1, replace=F)]
        lista_medoides_disponibles<-setdiff(lista_medoides_disponibles,nuevo_medoide)
        
        #puedo llegar a tener todos medoides explorados y quedarme sin medoide para asignar como nuevo medoide
        #por eso solo aplico cuando haya al menos un medoide disponible en lista_medoides_disponibles
        if(length(lista_medoides_disponibles)>0){
          
          #creo la solucion vecina S_prima numero i, esto es la misma que S pero con un gen alterado
          solucion_S_prima<-solucion_S
          solucion_S_prima[1,(posicion_alterar+1)]=nuevo_medoide
          
          #Verificar si la S prima produce singletons, si es asi no hago nada y paso a la sigueinte s_prima
          tablagroupsS_prima<-generate.groups(nrow(solucion_S_prima), as.matrix(solucion_S_prima[,2:(num_k+1)]), dmatrix1, dmatrix2)
          #quitar las soluciones que producen singletons. Actualiza con ello tabla_groups y poblacion
          arreglar<-singletons.delete(tablagroupsS_prima, solucion_S_prima, num_k)
          
          #preguntar si hay al menos una fila en cada lista.
          #
          #Si hay entonces seguir el proceso
          #sinO es porque solucion_s_prima produce singletons. No hago nada y paso a la sigueinte s-prima
          if(length(arreglar$groups)>0){
            
            
            #Evaluo cada solucion vecina
            poblacion_unida<-rbind(solucion_S_prima, poblacion_A)
            #Dejar soluciones sin repeticion
            poblacion_unida<-poblacion_unida[!duplicated(poblacion_unida[,2:(num_k+1)]),]
            
            #calcula dominancias de poblacion archivo y vecina S' (poblacion_N[i,])
            table.groups.PU<-generate.groups(nrow(poblacion_unida), as.matrix(poblacion_unida[,2:(num_k+1)]), dmatrix1, dmatrix2)
            #quitar las soluciones que producen singletons. Actualiza con ello tabla_groups y poblacion
            arreglar<-singletons.delete(table.groups.PU, poblacion_unida, num_k)
            poblacion_unida=arreglar$poblacion
            table.groups.PU=arreglar$groups
            
            
            
            dominanciasPU<-calculate.ranking.crowding(nrow(poblacion_unida), as.matrix(poblacion_unida[,1:(num_k+1)]), table.groups.PU, number.objectives, TRUE, dmatrix1, dmatrix2, num_k)
            dominanciasPU<-as.data.frame(dominanciasPU)
            
            #Une las columnas comunes (medoides en este caso) de dominancias (medoides y jerarquia) y s_prima
            calidad_solucion_S_prima<-merge(dominanciasPU, solucion_S_prima)
            
            #Si ninguna solucion del archivo dominada a la solucion vecina, o sea si solucion vecina es rnkIndex=1
            #Entonces agrego la solucion al archivo y elimino todas las soluciones que puedan estar dominada por ella
            #Si no es mejor o sea rnkIndex>=2 entonces ignoro esas s_prima y paso a la siguiente s_prima
            ###if(calidad_solucion_S_prima$rnkIndex=="1"){
            if(is.na(calidad_solucion_S_prima$crowding)) calidad_solucion_S_prima$crowding=0
            if(is.na(dominanciasPU$crowding[1])) dominanciasPU$crowding[1]=0
            
            if(calidad_solucion_S_prima$paretoranking=="1" && calidad_solucion_S_prima$crowding>=dominanciasPU$crowding[1]){
              
              
              #buscar solucion vecina S prima, y se marca como no explorada (0=FALSE)
              A=dominanciasPU[,2:(num_k+1)]
              B=solucion_S_prima[,2:(num_k+1)]
              
              if(is.vector(A)==TRUE) {
                A=as.matrix(A)
              }
              if(is.vector(B)==TRUE) {
                B=as.matrix(B)
              }
              
              fila_prima<-unname(which(apply(A, 1, function(x) all(x == B))))
              
              #buscar solucion vecina S prima, y se marca como no explorada (0=FALSE)
              #fila_prima<-unname(which(apply(dominanciasPU[,2:(num_k+1)], 1, function(x) all(x == solucion_S_prima[,2:(num_k+1)]))))
              
              
              #por si acaso este mas de una vez la solucion s prima, se marca como no explorada todas las veces que este
              for (x in 1: length(fila_prima)) {
                dominanciasPU[fila_prima[x],1]=0 #(0=FALSE)
              }
              
              solo_las_no_dominadas<-as.data.frame(subset(dominanciasPU, dominanciasPU$paretoranking=='1'))
              
              #actualiza la poblacion archivo  A dejando solo los no-dominados
              poblacion_A<-solo_las_no_dominadas[,1:(num_k+1)]
              
              ###solo_las_no_dominadas<-subset(dominanciasPU, rnkIndex=='1')
              #actualiza la poblacion archivo  A dejando solo los no-dominados
              ###poblacion_A<-solo_las_no_dominadas[,1:(num_k+1)]
              
            }
            
            
          }
          
        }
        
      }
      
      
      
      solucion_S<-as.data.frame(solucion_S)
      
      
      #buscar solucion S en poblacion A y marcar como explorada (1=TRUE)
      J=poblacion_A[,2:(num_k+1)]
      P=solucion_S[,2:(num_k+1)]
      if(is.vector(J)==TRUE) {
        J=as.matrix(J)
      }
      if(is.vector(P)==TRUE) {
        P=as.matrix(P)
      }
      
      fila_solucion_s<-unname(which(apply(J, 1, function(x) all(x == P))))
      
      
      if (length(fila_solucion_s)>0) {
        #por si acaso este mas de una vez la solucion s, marco como visitada todas las veces que este
        for (x in 1:length(fila_solucion_s)) {
          poblacion_A[fila_solucion_s[x],1]=1 #(1=TRUE)
        }
      }
      
      #Asigno a poblacion A0 Todas las soluciones S que estan en poblacion A, pero que no han sido visitadas
      poblacion_A<-as.data.frame(poblacion_A)
      
      # Quitar filas con NA
      poblacion_A0<-as.data.frame(subset(poblacion_A, poblacion_A[,1]=='0'))
      
      
      
    }
    
    #cat("\n  poblacion despues \n")
    #print(poblacion_A)
    
    #Dejar solo medoides en poblacion A
    poblacion_A<-poblacion_A[,2:(num_k+1)]
    colnames(poblacion_A)<-paste("V",1:(num_k),sep = "")
    
    
    
    #dejar sin cromosomas cromosomas repetidos
    poblacion_A<-poblacion_A[!duplicated(poblacion_A[,1:(num_k)]),]
    poblacion_A<-as.data.frame(poblacion_A)
    
    
    
    #calcula dominancias de poblacion archivo y vecina S' (poblacion_N[i,])
    ###tablagroupsA<-generate.groups(nrow(poblacion_A[1:pop_size,1:(num_k)]), as.matrix(poblacion_A[1:pop_size,1:(num_k)]), dmatrix1, dmatrix2, alfa)
    ###poblacion_A<-calculate.ranking.crowding(nrow(poblacion_A[1:pop_size,1:(num_k)]), as.matrix(poblacion_A[1:pop_size,1:(num_k)]), tablagroupsA, number.objectives, FALSE, matriz_exp, dmatrix1, dmatrix2, num_k,alfa, num_k)
    ###poblacion_A<-as.data.frame(poblacion_A)
    tablagroupsA<-generate.groups(nrow(poblacion_A[,1:(num_k)]), as.matrix(poblacion_A[,1:(num_k)]), dmatrix1, dmatrix2)
    poblacion_A<-calculate.ranking.crowding(nrow(poblacion_A[,1:(num_k)]), as.matrix(poblacion_A[,1:(num_k)]), tablagroupsA, number.objectives, FALSE, dmatrix1, dmatrix2, num_k)
    poblacion_A<-as.data.frame(poblacion_A)
    ###rownames(poblacion_A)<-c(1:nrow(poblacion_A))
    
    
    
    poblacion_jorge=subset(poblacion_A, poblacion_A$paretoranking==1)
    
    
    
    return(as.data.frame(poblacion_A[,]))
    
    
  }else{
    
    
    rownames(population_pareto)<-c(1:nrow(population_pareto))
    poblacion_A<-population_pareto[1,]
    return(poblacion_A)
    
  }
  
  
}




generate.path.relinking<-function(population_mejorar, num_k, dmatrix1, dmatrix2, number.objectives, cores){
  
  
  #population_mejorar=population.pareto
  #solo se generan combinaciones de a 2 si la poblacion a mejorar tiene mas de 2 soluciones
  if(nrow(population_mejorar)>1){
    
    if(nrow(population_mejorar)==2){
      
      #pares para PR(x,y)
      pares_normal=t(c(1,2))
      #pares para PR(y,x)
      pares_inverso=t(c(2,1))
      
    }else{
      #pares para PR(x,y)
      pares_normal=t(utils::combn(nrow(population_mejorar), 2))
      #pares para PR(y,x)
      pares_inverso=pares_normal[,c(2,1)]
    }
    
    
    
    #Univer ambos pares, o sea pares (x,y) con pares (y,x)
    pares<-rbind(pares_normal,pares_inverso)
    
    #almacenara las soluciones intermedias seleccionadas
    soluciones_path<-data.frame()
    iteraciones_path<-data.frame()
    #paralelismo
    message("empezar sin paralelismo")
    message(str_interp("Pares: ${nrow(pares)}"))
    for(i in 1:nrow(pares) ){
      message(i)
      if(i %% 100 == 0) {
        message(str_interp("${i}/${nrow(pares)}"))
      }
      
      invisible(rbind) ## Para solucionar problema de consumo excesivo de memoria . rbind is referenced and thus exported to workers
      # for (i in 1:nrow(pares)) {
      #paralelismo
      soluciones_path<-data.frame()
      
      s_initial=population_mejorar[pares[i,1],1:num_k]
      s_guide=population_mejorar[pares[i,2],1:num_k]
      
      
      #computar el numero de elementos en s_initial que no estan en s_guide
      #computar el numero de elementos en s_guide que no estan en s_initial
      #con lo anterior conocer? la diferencia simetrica entre s_initial y s_guide
      #es decir, los moves que faltan para llevar s_initial a s_guide
      
      while((length(setdiff(s_initial,s_guide)))!=0){
        
        #analizar todos los movimientos m, para llegar de s_initial a s_guide
        #eliminando un elemento de s_initial que no este en s_guide, y reemplazando
        #por un elemento de s_guide que no este en s_initial
        
        elementos_eliminar<-colnames(setdiff(s_initial,s_guide))
        elementos_agregar<-colnames(setdiff(s_guide,s_initial))
        
        soluciones_intermedias<-data.frame()
        intermedia=1
        
        for (mov in 1:(length(setdiff(s_initial,s_guide)))) {
          
          soluciones_intermedias<-rbind(soluciones_intermedias,s_initial)
          soluciones_intermedias[intermedia,c(elementos_eliminar[mov])]=s_guide[1,c(elementos_agregar[mov])]
          intermedia=intermedia+1
        }
        
        #calcular paretos y ver cual de las soluciones intermedias hasta el momento es la mejor
        #esta es seleccionada como la nueva s_initial y se almacena como una de las soluciones path
        #a ser consideradas.
        
        #si solo hay una intermedia significa que es la s_guide, por tanto no se agrega
        #a las soluciones path solo se pone como s_initial para que salga del while
        if(nrow(soluciones_intermedias)>1){
          
          tablagroupsIntermedias<-generate.groups(nrow(soluciones_intermedias), as.matrix(soluciones_intermedias[,1:(num_k)]), dmatrix1, dmatrix2)
          
          #quitar las soluciones que producen singletons. Actualiza con ello tabla_groups y poblacion
          arreglar<-singletons.delete(tablagroupsIntermedias, soluciones_intermedias, num_k)
          soluciones_intermedias=arreglar$poblacion
          tablagroupsIntermedias=arreglar$groups
          
          if(length(arreglar$groups)>0){
            dominanciasIntermedias<-calculate.ranking.crowding(nrow(soluciones_intermedias), as.matrix(soluciones_intermedias[,1:(num_k)]), tablagroupsIntermedias, number.objectives, FALSE, dmatrix1, dmatrix2, num_k)
            dominanciasIntermedias<-as.data.frame(dominanciasIntermedias)
            
            soluciones_path<-rbind(soluciones_path, dominanciasIntermedias[1, 1:num_k])
            s_initial=dominanciasIntermedias[1,1:num_k]
            
          }else{
            #como la intermedias tiene soluciones que producen singleton, no considero ninguna y
            #para que salga del while asigno s_guide como s_initial
            s_initial=s_guide
          }
          
          
        }else{
          
          #ver si la unica solucion en soluciones_intermedia es singleton.
          #de ser afirmativo asigno s_guide como s_initial para que salga del while
          tablagroupsIntermedias<-generate.groups(nrow(soluciones_intermedias), as.matrix(soluciones_intermedias[,1:(num_k)]), dmatrix1, dmatrix2)
          #quitar las soluciones que producen singletons. Actualiza con ello tabla_groups y poblacion
          arreglar<-singletons.delete(tablagroupsIntermedias, soluciones_intermedias, num_k)
          soluciones_intermedias=arreglar$poblacion
          tablagroupsIntermedias=arreglar$groups
          
          if(length(arreglar$groups)>0){
            s_initial=soluciones_intermedias[1,1:num_k]
          }else{
            s_initial=s_guide
          }
          
        }
        
        
      }
      
      
      #almacena para cada i las soluciones path
      iteraciones_path=rbind(iteraciones_path, soluciones_path)
      
    }
    message("BL sin paralelismo finalizada")
    
    
    #Unir las soluciones path, con las soluciones elite (o sea las pareto)
    #Determinar cuales son las que definitivamente quedan mejoradas, es decir las no-dominadas de todas estas
    
    #a veces solo hay dos soluciones pareto que tienen solo un gen diferente, eso
    #ocasiona que no se agregre nada a la variable soluciones_path, por ejemplo
    #87 123 34 12
    #87 123 31 12
    #por ello si hago
    #rbind da error, as?,  solo se hace el rbind cuando no pasa esto, en otro caso el merge es solo la
    #poblacion a mejorar
    
    #cat("soluciones en path", nrow(iteraciones_path), "\n")
    #print(iteraciones_path)
    
    if (nrow(iteraciones_path)>0) {
      #paralelismo
      
      #soluciones_path[, 1:num_k]=as.matrix(iteraciones_path)
      soluciones_path=as.matrix(iteraciones_path)
      #print(soluciones_path)
      
      soluciones_merge<-rbind(population_mejorar[,1:num_k], soluciones_path[, 1:num_k])
      
    }else{
      soluciones_merge<-population_mejorar[,1:num_k]
      
    }
    
    
    tablagroupsMerge<-generate.groups(nrow(soluciones_merge), as.matrix(soluciones_merge[,1:(num_k)]), dmatrix1, dmatrix2)
    
    #Quitar soluciones con singletons o con cluster sin genes.
    arreglar<-singletons.delete(tablagroupsMerge, as.matrix(soluciones_merge[,1:(num_k)]), num_k)
    
    #cat("arreglar", nrow(arreglar$poblacion),"\n")
    
    soluciones_merge=arreglar$poblacion
    tablagroupsMerge=arreglar$groups
    
    if(length(arreglar$groups)>0){
      
      dominanciasMerge<-calculate.ranking.crowding(nrow(soluciones_merge), as.matrix(soluciones_merge[,1:(num_k)]), tablagroupsMerge, number.objectives, FALSE, dmatrix1, dmatrix2, num_k)
      dominanciasMerge<-as.data.frame(dominanciasMerge)
      
      
      #print(dominanciasMerge)
      
      soluciones_no_dominadas<-subset(dominanciasMerge, dominanciasMerge$paretoranking=='1')
      
      soluciones_no_dominadas<-soluciones_no_dominadas[!duplicated(soluciones_no_dominadas[,1:(num_k)]),]
      #cat("\n")
      #cat( " no-dominadas ", nrow(soluciones_no_dominadas), "\n")
      rownames(soluciones_no_dominadas)<-c(1:nrow(soluciones_no_dominadas))
      #print(soluciones_no_dominadas)
      
      
    }else{
      soluciones_no_dominadas=population_mejorar
    }
    
    
    return(soluciones_no_dominadas)
    
    
    
  }else{
    
    warning("It is not possible to apply Path-Relinking (Only a solution) ...")
    rownames(population_mejorar)<-c(1:nrow(population_mejorar))
    return(population_mejorar)
    
    
  }
  
  
  
  
}


generate.results<-function(num_k, dmatrix1, dmatrix2, pop_size, generation, rat_cross, rat_muta, population){
  
  
  ############ obtener Pareto de la ultima generacion
  d<-as.data.frame(population[1:pop_size,])
  D = d[order(d$paretoranking, -d$crowding), ]
  #las que son Pareto
  paretos<-subset(D, D$paretoranking=='1')
  nondom<-paretos
  
  #Forma groups de cada solucion Pareto
  population.PARETO<-as.matrix(nondom[,1:num_k])
  table.groups.PARETO<-generate.groups(nrow(population.PARETO), population.PARETO, dmatrix1, dmatrix2)
  
  #Guardar soluciones pareto
  pareto_solutions<-as.data.frame(table.groups.PARETO)
  colnames(pareto_solutions)<-1:length(table.groups.PARETO)
  
  vectors.partition.list=list()
  for (i in 1:ncol(pareto_solutions)) {
    vector.partition=unlist(pareto_solutions[, i])
    names(vector.partition)<-rownames(pareto_solutions)
    vectors.partition.list[[i]]=vector.partition
  }
  
  
  
  return(list("population"=nondom,
              "matrix.solutions"=pareto_solutions,
              "clustering"=vectors.partition.list))
  
}







#' Perform the Multi-Objective Clustering Algorithm Guided by a-Priori Biological Knowledge (MOC-GaPBK)
#'
#' @description This function receives two distance matrices and it performs the MOC-GaPBK.
#' @param dmatrix1 A distance matrix. It should have the same dimensions that dmatrix2. It is mandatory.
#' @param dmatrix2 A distance matrix. It should have the same dimensions that dmatrix1. It is mandatory.
#' @param num_k  The number k of groups represented by medoids in each individual. It is mandatory.
#' @param generation Number of generations to be performed by MOC-GaPBK. By default 50.
#' @param pop_size Size of population. By default 10.
#' @param rat_cross Probability of crossover. By default 0.80.
#' @param rat_muta Probability of mutation. By default 0.01.
#' @param tour_size Size of tournament. By default 2.
#' @param neighborhood Percentage of neighborhood. A real value between 0 and 1. It is computed as neighborhood*pop_size to determine the size of neighborhood. By default 0.10.
#' @param local_search A boolean value indicating whether the local searches procedures (PR and PLS) are computed. By default \emph{FALSE}.
#' @param cores Number of cores to be used to compute the local searches procedures. By default 2.
#' @return
#' \item{population}{The population of medoids including the objective functions values and order by Pareto ranking and crowding distance values.}
#' \item{matrix.solutions}{A matrix with results of clustering. Each column represents a clustering solution available in Pareto front.}
#' \item{clustering}{A list containing named vectors of integers from 1:k representing the cluster to which each object is assigned.}
#' @author  Jorge Parraga-Alava, Marcio Dorn, Mario Inostroza-Ponta
#' @keywords GaPBK
#' @keywords NSGA-II
#' @keywords Pareto ranking
#' @keywords Crowding distance
#' @keywords Xie-Beni validity index
#' @details MOC-GaPBK is a method proposes by Parraga-Alava, J. et. al. 2018. It carries out the discovery of clusters using NSGA-II algorithm along with Path-Relinking (PR) and Pareto Local Search (PLS) as intensification and diversification strategies, respectively. The algorithm uses as objective functions two versions of the Xie-Beni validity index, i.e., a version for each distance matrix (dmatrix1, dmatrix2). More details about this compute can be found in: <https://doi.org/10.1186/s13040-018-0178-4>. MOC-GaPBK yield a set of the best clustering solutions from a multi-objective point of views.
#' @export
#' @import foreach
#' @import doParallel
#' @import doMPI
#' @import doSNOW
#' @import parallel
#'
#' @examples
#'
#' ##Generates a data matrix of dimension 50X20
#'
#' library("amap")
#' library("moc.gapbk")
#'
#' x <- matrix(runif(100*20, min = -5, max = 10), nrow=50, ncol=20)
#'
#' ##Compute two distance matrices
#'
#' dmatrix1<- as.matrix(amap::Dist(x, method = "euclidean"))
#' dmatrix2<- as.matrix(amap::Dist(x, method = "correlation"))
#'
#' ##Performs MOC-GaPBK with 5 cluster
#'
#' example<-moc.gabk(dmatrix1, dmatrix2, 5)
#'
#' example$population
#' example$matrix.solutions
#' example$clustering
#'
#' @references
#' J. Parraga-Alava, M. Dorn, M. Inostroza-Ponta (2018). \emph{A multi-objective gene clustering algorithm guided by apriori biological knowledge with intensification and diversification strategies}. BioData Mining. 11(1) 1-16.
#'
#' K. Deb, A. Pratap, S. Agarwal, T. Meyarivan (2002). \emph{A fast and elitist multiobjective genetic algorithm: NSGA-II}. IEEE Transactions on Evolutionary Computation, 6(2) 182-197.
#'
#' F. Glover (1997). \emph{Tabu Search and Adaptive Memory Programming - Advances, Applications and Challenges}. Interfaces in Computer Science and Operations Research: Advances in Metaheuristics, Optimization, and Stochastic Modeling Technologies. 1-75.
#'
#' J. Dubois-Lacoste, M. Lopez-Ibanez, Stutzle, T. (2015). \emph{Anytime Pareto local search}. European Journal of Operational Research, 243(2) 369-385.

moc.gapbk<-function(dmatrix1, dmatrix2, num_k,
                   generation=50, pop_size=10, rat_cross=0.80, rat_muta=0.01, tour_size=2,
                   neighborhood=0.10, local_search=FALSE, cores=2){
  
  
  if(base::nrow(dmatrix1)==base::ncol(dmatrix1) && base::nrow(dmatrix2)==base::ncol(dmatrix2)){
    if(base::nrow(dmatrix1)==base::nrow(dmatrix2) && base::ncol(dmatrix1)==base::ncol(dmatrix2)){
      if(exists('num_k')==TRUE){
        if(num_k>0){
          
          if(neighborhood>=0 & neighborhood<=1){
            
            #Set parameters
            num_objects=nrow(dmatrix1)
            number.objectives=2
            
            #Initialization
            population.P=generate.initial.population(nrow(dmatrix1), num_k, pop_size)
            
            g=1
            
            while (g<=generation) {
              message(str_interp("Generation ${g}/${generation}"))
              message(Sys.time())
              # Population P
              verify.singletons.P<- singletons.repair(population.P, dmatrix1, dmatrix2, num_objects)
              table.groups.P<- verify.singletons.P[[1]] #1 corresponde al primer argumento devuelto y corresponde a groups formados
              population.P<- as.matrix(verify.singletons.P[[2]]) #1 corresponde al segundo argumento ... es la poblacion sin singletons
              population.P<- calculate.ranking.crowding(pop_size, population.P, table.groups.P, number.objectives, FALSE, dmatrix1, dmatrix2, num_k)
              #selection
              mating.pool <- nsga2R::tournamentSelection(population.P, pop_size, tour_size)
              population.Q <- t(sapply(1:pop_size, function(u) array(rep(0,num_k))))
              #crossover
              crossover <-generate.crossover.k.points(pop_size,num_k, mating.pool, population.Q, rat_cross)
              population.Q<-verify.feasibility(crossover,population.Q, num_objects)
              #mutation
              mutation <-generate.mutation.random.controller(pop_size, num_k, population.Q, rat_muta, num_objects)
              population.Q<-verify.feasibility(mutation, population.Q, num_objects)
              # Population Q
              verify.singletons.Q<-singletons.repair(population.Q, dmatrix1, dmatrix2, num_objects)
              table.groups.Q<-verify.singletons.Q[[1]]
              population.Q<- as.matrix(verify.singletons.Q[[2]])
              population.Q<-calculate.ranking.crowding(pop_size, population.Q, table.groups.Q, number.objectives, FALSE, dmatrix1, dmatrix2, num_k)
              # Population R
              #como P y Q han sido corregidas para no tener singletons, ya no debo verificar eso en poblacion R, ya que es la union de ambas
              population.R <- rbind(population.P, population.Q)
              rownames(population.R)<-1:nrow(population.R)
              
              population.R<-population.R[, -((num_k+1):(num_k+number.objectives+2))] #deja solo cromosomas, sea quita objetivos, Jerarq y crowding
              
              #recalcula Jerarquia Pareto y crowding en poblacion R
              table.groups.R<-generate.groups(nrow(population.R), population.R, dmatrix1, dmatrix2)
              population.R<-calculate.ranking.crowding(pop_size*2, population.R, table.groups.R, number.objectives, FALSE, dmatrix1, dmatrix2, num_k)
              population.R<-as.data.frame(population.R)
              #Solo usa soluciones Pareto para las busqueda locales y para rellenar la poblacion.
              population.pareto<-subset(population.R, population.R$paretoranking=='1')
              population.pareto<-population.pareto[!duplicated(population.pareto[,1:(num_k)]),]
              population.pareto<-population.R[!duplicated(population.R[,1:(num_k)]),]
              #Quitar soluciones con singletons o con cluster sin genes.
              table.groups.Pareto<-generate.groups(nrow(population.pareto), as.matrix(population.pareto[,1:num_k]), dmatrix1, dmatrix2)
              arreglar<-singletons.delete(table.groups.Pareto, population.pareto, num_k)
              population.pareto=arreglar[[2]]
              population.R<-population.pareto
              
              if(local_search==TRUE){
                message("Comenzando búsqueda local (path-relinking)")
                #Path relinking
                population.R<-generate.path.relinking(population.pareto, num_k, dmatrix1, dmatrix2, number.objectives, cores)#, lista_funciones, lista_paquetes, cores) #devuelve poblacion R pero solo  Paretos sin repeticiones
                #Pareto Local Search
                message("Comenzando búsqueda local (pareto local search)")
                population.R<-generate.pareto.local.search(population.R, neighborhood, num_k, num_objects, pop_size,
                                                           dmatrix1, dmatrix2, number.objectives) #Recibe solo Paretos (devueltos en population.R por PR) con Rnk y crowding
                message("Búsqueda local finalizada")
              }
              
              # Population P+1
              #population.pareto=subset(population.R, population.R$paretoranking==1)
              if(nrow(population.R)< pop_size) {
                population.P<- as.matrix(population.R[,1:num_k])
              }else{
                population.P<- as.matrix(population.R[1:pop_size,1:num_k])
              }
              
              #Rellena poblacion con Pr
              if(nrow(population.P)< pop_size){
                population_random<-as.data.frame(t(sapply(1:pop_size, function(u) array(sample(1:num_objects, num_k, replace=F)))))
                #por si acaso se prodyjeron singletons en la poblacion ramdom
                reparar<-singletons.repair(as.matrix(population_random), dmatrix1, dmatrix2, num_objects)
                population_random=reparar[[2]] #2 porq el primero son groups formados
                #dejar sin cromosomas cromosomas repetidos
                if(nrow(population.P)>1){
                  population.P<-population.P[!duplicated(population.P[,1:(num_k)]),]
                  population.P<-as.data.frame(population.P)
                }
                population.P<-rbind(population.P[,1:(num_k)], population_random[(nrow(population.P)+1):pop_size,])
                rownames(population.P)<-c(1:nrow(population.P))
                population.P<-as.matrix(population.P)
              }
              
              
              g=g+1
              
            }
            
            message("Generando resultados...")
            results=generate.results(num_k, dmatrix1, dmatrix2, pop_size, generation, rat_cross, rat_muta, population.R)
            
            return(results)
            
            
          }else{
            stop("'neighborhood' should be between 0 and 1")
          }
          
        }else{
          stop("'k' should be more than 1")
        }
      }else{
        stop("'k' should be indicated")
      }
    }else{
      stop('Both matrices or data.frames should have equal number of rows and columns')
    }
  }else{
    stop("'dmatrix1' or 'dmatrix2' must have the same number of rows and columns")
  }
  
  
  
  
}