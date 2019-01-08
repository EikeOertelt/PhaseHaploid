###########################################################################################################

# Functions

###########################################################################################################

# main
phase = function(data, ExprBeginAtCol, DistanceCol, MinDensity=0.2, DistributionVaryBy=0.3, IncludeC=T, RemoveArtifacts=T){

  distribution = rcountAB(data)
  density = get_density(distribution, length(data[1,ExprBeginAtCol : length(data[1,])]))
  selection = data[c(1,which(((density > MinDensity) & (distribution[,1]
                                                        < (1 + DistributionVaryBy)) & (distribution[,1]  > (1 - DistributionVaryBy))))),]
  result = phasing(selection, ExprBeginAtCol)

  if(RemoveArtifacts == TRUE) {result = remove_artifacts(result, ExprBeginAtCol)}
  if(IncludeC == TRUE) {result = includeC(result, ExprBeginAtCol, DistanceCol)}

  result = cbind(result, cbind(get_distance(result),countrecomb(result, ExprBeginAtCol)))

  utils::write.table(result, file  = "result.csv", row.names=FALSE, col.names=FALSE, sep=";")

  return("Job complete, file saved")
}

phase_two = function(data1, data2, ExprBeginAtCol, DistanceCol, MinDensity=0.2, DistributionVaryBy=0.3, IncludeC=T, RemoveArtifacts=T){

  distribution1 = rcountAB(data1)
  density1 = get_density(distribution1, length(data1[1,ExprBeginAtCol : length(data1[1,])]))
  selection1 = data1[c(1,which(((density1 > MinDensity) & (distribution1[,1]
                                                           < (1 + DistributionVaryBy)) & (distribution1[,1]  > (1 - DistributionVaryBy))))),]

  distribution2 = rcountAB(data2)
  density2 = get_density(distribution2, length(data2[1,ExprBeginAtCol : length(data2[1,])]))
  selection2 = data2[c(1,which(((density2 > MinDensity) & (distribution2[,1]
                                                           < (1 + DistributionVaryBy)) & (distribution2[,1]  > (1 - DistributionVaryBy))))),]

  result = phasing_and_merge(selection1, selection2, ExprBeginAtCol)


  if(RemoveArtifacts == TRUE) {result = remove_artifacts(result,ExprBeginAtCol)}
  if(IncludeC == TRUE) {result = includeC(result, ExprBeginAtCol,DistanceCol)}

  result = cbind(result, cbind(get_distance(result),countrecomb(result, ExprBeginAtCol)))

  utils::write.table(result, file  = "result.csv", row.names=FALSE, col.names=FALSE, sep=";")

  return("Job completet, file saved")
}


# preprocessing
rcountAB = function(dataset) {
  # this function counts the occurances of the 2 phases for each row and returns a 3 column matrix
  # first column: the A to B ratio
  # second column: number of As
  #thrid column: number of Bs


  # a matrix that will later carry the result
  # three coloumns and a rows depending on the length of the dataset
  erg = matrix(0,length(as.matrix(dataset[,1])),3)

  # we need to transform the dataframe into a matrix, because dataframes tend to behave strange
  dataset = as.matrix(dataset)

  for(i in 1:length(dataset[,1])) {
    # count the As and Bs
    nrA  = length(which(grepl("1", dataset[i,])))
    nrB  = length(which(grepl("0", dataset[i,])))

    # if there are no Bs we need to manually set the value of this row to catch an error
    if(nrB == 0) {
      erg[i,1] = 0
      erg[i,2] = nrA
      erg[i,3] = nrB
      next
    }

    erg[i,1] = nrA/nrB
    erg[i,2] = nrA
    erg[i,3] = nrB
  }
  return(erg)
}

get_density = function (dataset, nrExpr) {
  # calculates the density of markers in a row
  # therefor it needs the result of the rcountAB function

  erg = matrix(0,length(dataset[,1]), 1)
  end = length(dataset[,1])
  lastc = length(dataset[1,])
  for (i in  1:end) {
    erg[i]= (as.numeric(dataset[i,lastc]) + as.numeric(dataset[i, lastc-1])) / nrExpr
  }
  return (erg)
}

get_distance = function(dataset) {
  # calculates the dictances of each position to the previous position

  dataset = as.matrix(dataset)
  erg = matrix(0,length(dataset[,1]),1)
  erg[1] = "DISTANCE"

  NrLGs = length(as.character(unique(dataset[,1]))) #number of LGs

  for(k in 1:NrLGs) {
    vec = which(dataset[,1] == as.character(unique(dataset[,1])[k]))

    if(length(vec)>1){
      for(i in vec[-1]) {
        erg[i,1] = as.numeric(dataset[i,2]) - as.numeric(dataset[i-1,2])
      }
    }
  }
  return(erg)
}


# phasing
initialise = function(datarow, colstart, A, B) {
  # changes all cells in a given row to the phase it belongs

  end = length(datarow)
  datarow = as.matrix(datarow)
  for ( i in colstart:end) {
    if(grepl(A, as.character(datarow[i]))) {
      datarow[i] = "A"
      next
    }
    if(grepl(B, as.character(datarow[i]))) {
      datarow[i] = "B"
      next
    }
    if(grepl("[.]", as.character(datarow[i]))) {
      datarow[i] = "-"
      next
    }
  }
  return (datarow)
}

get_phase = function(dataset, row, AB) {
  # tries to infer the phase for a row
  # to do this it compares a line with the allready phased previous line

  eq = 0 # the number of witnesses that the second line has the same AB as the first
  neq = 0 # the number of witnesses that the second line has NOT the same AB as the first

  # compares the cells of the previous row with the ones from the actual row
  # depending on which expression we find either eq or neq is increased
  for(i in which(grepl(AB[1], as.matrix(dataset[row-1, ])))) {
    if(grepl(AB[1], dataset[row, i])) {
      eq = eq + 1
    }
    if(grepl(AB[2], dataset[row, i])) {
      neq = neq + 1
    }
  }

  if(eq > neq) {
    return(c(AB[1],AB[2]))
  }
  return(c(AB[2],AB[1]))
}

get_phase_reversed = function(dataset, row, AB) {
  # does the same as the above function, just in the opposite direction

  eq = 0
  neq = 0
  for(i in which(grepl(AB[1], as.matrix(dataset[row+1, ])))) {
    if(grepl(AB[1], dataset[row, i])) {
      eq = eq + 1
    }
    if(grepl(AB[2], dataset[row, i])) {
      neq = neq + 1
    }
  }

  if(eq > neq) {
    return(c(AB[1],AB[2]))
  }
  return(c(AB[2],AB[1]))
}

find_first_shared_loci = function(dataset1, dataset2, vec1, vec2) {
  # searches for the first shared loci in the two datasets
  # gets two vectors containing all rows for a LG and the datasets themself
  # it just checks whether two positions are the same in this LG
  erg = intersect(dataset1[vec1,2], dataset2[vec2,2])

  #if there is no shared loci, return 0
  if(length(erg)==0) {return(0)}

  # return 0 if there s no shared locus
  return(max(erg[1],0))
}

phasing_and_merge = function(dataset1, dataset2, colstart) {
  # main function that applies the other function on the datasets
  # merges both dataset1s at the end

  dataset1 = as.matrix(dataset1)
  dataset2 = as.matrix(dataset2)

  NrLGs = getLGcount(dataset1) # number of LGs (equal for both datasets)
  # it is expected, that the first row is a heading

  # initialisation of the two intermediate results
  # maybe not the best performance, because we copy two huge datasets, but the easiest solution
  zwerg1 = dataset1
  zwerg2 = dataset2


  # this loop iterates all the LG (excluding the first, which is supposed to be a heading)
  for(k in 2:NrLGs){
    # first we extract the rows belonging to a LG
    vec1 = getVecOfLG(zwerg1, k) # vector of the given LG in dataset1
    vec2 = getVecOfLG(zwerg2, k) # vector of the given LG in dataset2

    # now we look for the first shared locus
    # this will later be our anchor point in the phasing
    start1 = which(zwerg1[,2] == find_first_shared_loci(zwerg1, zwerg2, vec1, vec2)
                   & (zwerg1[vec1[1],1] == unique(zwerg2[,1])[k]))
    if(length(start1)==0) {start1 = vec1[1]}

    start2 = which(zwerg2[,2] == find_first_shared_loci(zwerg1, zwerg2, vec1, vec2)
                   & (zwerg2[vec2[1],1] == unique(zwerg2[,1])[k]))
    if(length(start2)==0) {start2 = vec2[1]}

    # we started with A as 1/1 and B as 0/0
    # this variable belongs to the initialisation, it tells which value will be translated to A or rather B
    AB = c("1", "0")

    #initialising each LG at the shared loci
    zwerg1[start1, ] = initialise(dataset1[start1, ], colstart, AB[1], AB[2])
    zwerg2[start2, ] = initialise(dataset2[start2, ], colstart, AB[1], AB[2])

    #if the first loci is the shared one, we cant look up
    if(vec1[1] != start1) {
      AB = c("1", "0")
      for(i in (start1-1) : vec1[1]){
        AB = get_phase_reversed(dataset1, i, AB)
        zwerg1[i, ] = initialise(zwerg1[i, ], colstart, AB[1], AB[2])
      }
    }

    #if vec1 has just 1 entry, we allready dealt with it
    if(length(vec1) > 1) {
      end = vec1[length(vec1)]
      AB = c("1", "0")
      for(i in (start1+1) : end){
        AB = get_phase(dataset1, i, AB)
        zwerg1[i, ] = initialise(zwerg1[i, ], colstart, AB[1], AB[2])
      }
    }


    #and all the same for the second dataset
    if(vec2[1] != start2) {
      AB = c("1", "0")
      for(i in (start2-1) : vec2[1]){
        AB = get_phase_reversed(dataset2, i, AB)
        zwerg2[i, ] = initialise(zwerg2[i, ], colstart, AB[1], AB[2])
      }
    }


    if(length(vec2) > 1) {
      end = vec2[length(vec2)]
      AB = c("1", "0")
      for(i in (start2+1) : end){
        AB = get_phase(dataset2, i, AB)
        zwerg2[i, ] = initialise(zwerg2[i, ], colstart, AB[1], AB[2])
      }
    }
  }

  # now the two datasets get finally merged
  erg = merge_two_datasets(zwerg1, zwerg2, colstart)
  return(erg)
}

phasing = function(dataset, colstart) {
  # main function that applies the other function on the datasets

  dataset = as.matrix(dataset)

  NrLGs = getLGcount(dataset) # number of LGs
  # it is expected, that the first row is a heading

  # initialisation of the intermediate result
  zwerg = dataset


  # this loop iterates all the LG (excluding the first, which is supposed to be a heading)
  for(k in 2:NrLGs){
    # first we extract the rows belonging to a LG
    vec = getVecOfLG(zwerg, k) # vector of the given LG in dataset

    # this will later be our anchor point in the phasing
    start = vec[1]

    # we started with A as 1/1 and B as 0/0
    # this variable belongs to the initialisation, it tells which value will be translated to A or rather B
    AB = c("1", "0")

    #initialising each LG at the shared loci
    zwerg[start, ] = initialise(dataset[start, ], colstart, AB[1], AB[2])

    # the calculation of 'end' boosts the performance of the loop
    # otherwise it would have to calculate the end at each step of the iteration again
    end = vec[length(vec)]

    AB = c("1", "0")
    if(start < end) {
      for(i in (start+1) : end){
        AB = get_phase(dataset, i, AB)
        zwerg[i, ] = initialise(zwerg[i, ], colstart, AB[1], AB[2])
      }
    }
  }

  return(zwerg)
}

merge_two_datasets = function(dataset1, dataset2, ExprBeginAtCol){
  # merges two datasets into one
  # they need to have the same format (start of expressiontable, etc. )

  NrLGs = length(as.character(unique(dataset1[,1]))) # number of LGs (equal for both datasets)

  # the first row of the result is the merging of heading from the first and all the individuals
  # from the second dataset
  erg = matrix("",length(dataset1[,1]) + length(dataset2[,1]), length(dataset1[1, ]) + length(dataset2[1, ExprBeginAtCol:length(dataset2[1,])]))
  erg[1, ] = c(dataset1[1, ], dataset2[1,ExprBeginAtCol:length(dataset2[1,])])

  #global row-counter, since we are using a while-loop
  i = 1

  # this loop iterates all the LG (excluding the first, which is supposed to be a heading)
  for(k in 2:NrLGs){

    # first we extract the rows belonging to a LG
    vec1 = which(dataset1[,1] == as.character(unique(dataset1[,1])[k])) # vector of the given LG in dataset1
    vec2 = which(dataset2[,1] == as.character(unique(dataset2[,1])[k])) # vector of the given LG in dataset2

    current1 = 1 # the current row for the first dataset
    current2 = 1 # the current row for the second dataset


    #calculating this befor the loop boosts the performance
    lengthvec1 = length(vec1)
    lengthvec2 = length(vec2)


    # loop that works till one of the datasets is at the end of its LG
    while((current1 <= lengthvec1) & (current2 <= lengthvec2)){
      i = i+1

      # first case: position in dataset1 < pos in dataset2
      if(as.numeric(dataset1[vec1[current1],2]) < as.numeric(dataset2[vec2[current2],2])) {

        # merge the two data frames and bind them to the result
        erg[i, ] = c(dataset1[vec1[current1], ], matrix("-", 1, length(dataset2[1,ExprBeginAtCol:length(dataset2[1,])])))
        current1 = current1 + 1
        next
      }


      # second case: position in dataset1 > pos in dataset2
      if(as.numeric(dataset1[vec1[current1],2]) > as.numeric(dataset2[vec2[current2],2])) {

        # merge the two data frames and bind them to the result
        erg[i, ] = c(dataset2[vec2[current2],1:(ExprBeginAtCol-1)], matrix("-", 1, length(dataset1[1,ExprBeginAtCol:length(dataset1[1,])])), dataset2[vec2[current2], ExprBeginAtCol :length(dataset2[1,])])
        current2 = current2 + 1
        next
      }


      # third case: position in dataset1 = pos in dataset2
      if(as.numeric(dataset1[vec1[current1],2]) == as.numeric(dataset2[vec2[current2],2])) {

        # merge the two data frames and bind them to the result
        erg[i, ] = c(dataset1[vec1[current1], ], dataset2[vec2[current2], ExprBeginAtCol :length(dataset2[1,])])
        current1 = current1 + 1
        current2 = current2 + 1
        next
      }
    }

    # bind the rest of dataset1 to the result
    if(current1 < lengthvec1) {
      for(j in current1 : lengthvec1){
        i = i + 1
        erg[i, ] = c(dataset1[vec1[j], ], matrix("-", 1, length(dataset2[1,ExprBeginAtCol:length(dataset2[1,])])))
      }
    }

    #bind the rest of dataset2 to the result
    if(current2 < lengthvec2) {
      for(j in current2 : lengthvec2){
        i = i + 1
        erg[i, ] = c(dataset2[vec2[j],1:(ExprBeginAtCol-1)], matrix("-", 1, length(dataset1[1,ExprBeginAtCol:length(dataset1[1,])])), dataset2[vec2[j], ExprBeginAtCol:length(dataset2[1,])])
      }
    }
  }
  return (erg)
}

remove_artifacts = function(dataset, ExprBeginAtCol){
  # this function removes the artifacts in a dataset, as artifact counts:
  # each row which has cells like 0/2, 1/2, 2/2, 0/3, 1/3, 2/3 or 3/3 in it

  # change the dataframe to a matrix
  dataset = as.matrix(dataset)
  len = length(dataset[1,])

  # result vector, saves all lines to delete
  erg = vector("integer", length(dataset[,1]))
  rowcounter = 1

  # iterate the dataset and save all lines with artifacts
  end = length(dataset[,1])
  for(i in  2:end){

    # now following: "best of copy paste" by eike oertelt
    # if we have at least one artifact in a row, the row will be saved and later deleted
    if(TRUE %in% grepl("0/2", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      rowcounter = rowcounter + 1
      next
    }

    if(TRUE %in% grepl("1/2", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      rowcounter = rowcounter + 1
      next
    }

    if(TRUE %in% grepl("2/2", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      rowcounter = rowcounter + 1
      next
    }

    if(TRUE %in% grepl("0/3", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      rowcounter = rowcounter + 1
      next
    }

    if(TRUE %in% grepl("1/3", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      rowcounter = rowcounter + 1
      next
    }

    if(TRUE %in% grepl("2/3", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      next
    }

    if(TRUE %in% grepl("3/3", dataset[i, ExprBeginAtCol:len])){
      erg[rowcounter] = i
      rowcounter = rowcounter + 1
      next
    }
  }
  #break if no artifacts are there
  if(sum(erg)==0){return(dataset)}

  dataset = dataset[-erg, ]
}


# misc
getLGcount = function(data) {
  # counts all the LG in a dataset
  erg = length(as.character(unique(data[,1])))
  return (erg)
}

getVecOfLG = function (data, k){
  #return a vector of all the loci belonging to a LG (k) in this Dataset
  erg = which(data[,1] == as.character(unique(data[,1])[k]))
  return(erg)
}

countrecomb = function(data, ExprBeginAtCol) {
  # counts the visual recombination events in a row

  # result: a vec of the numbers
  erg = matrix(0,length(data[,1]),1)
  erg[1] = "RECOMB.EVENTS"

  # for each LG
  end1 = getLGcount(data)
  for(k in 2 : end1){

    vec = getVecOfLG(data, k)

    if(length(vec) > 1) {
      # for all rows in this LG
      end2 = length(vec)
      for(i in 2:end2){

        #intermediate result counts the number of crossing overs
        zwerg = 0

        # for each cells in this row
        end3 = length(data[1,])
        for(j in ExprBeginAtCol : end3){

          # count recomb.events

          #if we have an A or a B in the current cell
          if((as.character(data[vec[i],j]) == "A") | (as.character(data[vec[i],j]) == "B")){

            #if the previous cell is "-" we have to check the first upstream cell that isnt a "-"

            bla = lookup(data, vec, i-1, j)
            if(bla != 0) {
              if(as.character(data[vec[i],j]) != (as.character(data[bla,j]))) {
                zwerg = zwerg + 1
              }
            }
          }
        }
        erg[vec[i]] = zwerg
      }
    }
  }
  return(erg)
}

lookup = function(data, vec, current, col) {
  # checks starting from the current position the upstream cells for their phase and returns the row
  # returns 0 if there is no A or B inside the index

  #current means the current position in the 'vec'-vector
  ABpos = 0

  for(i in (current):1){
    if((as.character(data[vec[i], col])=="A") | (as.character(data[vec[i], col])=="B" )){
      ABpos = vec[i]
      break
    }
  }

  return(ABpos)
}

lookdown = function(data, vec, current, col) {
  # checks starting from the current position the downstream cells for their phase and returns the row
  # returns 0 if there is no A or B inside the index

  ABpos = 0

  for(i in current:length(vec)){
    if((as.character(data[vec[i], col])=="B") | (as.character(data[vec[i], col])=="A" )){
      ABpos = vec[i]
      break
    }
  }

  return(ABpos)
}

includeC = function(data, ExprBeginAtCol, DistanceCol) {
  # we include the C by checking up- and downstrem (inside the column)

  # for each LG
  for(k in 2: (getLGcount(data))){

    vec = getVecOfLG(data, k)

    if(length(vec) > 1) {
      #for each row
      for(i in 1:length(vec)) {

        #for each cell
        end = length(data[1,])
        for(j in ExprBeginAtCol : end) {

          # check wheter the cell is a C
          if(grepl("C", as.character(data[vec[i],j]))) {

            # up and down give the row the first occurance of an A or B in the dataset
            up = lookup(data, vec, i, j)
            down = lookdown(data, vec, i, j)

            sumup = 0
            sumdown = 0

            # if up == 0, then there is no A or B and we dont need to deal with it
            if(up == 0) {
              sumup = Inf
            }
            else{
              sumup = as.numeric(data[vec[i], DistanceCol]) - as.numeric(data[up, DistanceCol])
            }

            # if down == 0, then there is no A or B and we dont need to deal with it
            if(down == 0){
              sumdown = Inf
            }
            else{
              sumdown = as.numeric(data[down, DistanceCol]) - as.numeric(data[vec[i], DistanceCol])
            }


            if(sumup > sumdown) {
              data[vec[i],j] = data[down,j]
            }
            if(sumup < sumdown) {
              data[vec[i],j] = data[up,j]
            }
          }
        }
      }
    }
  }

  return(data)
}
