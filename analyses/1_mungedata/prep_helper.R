

# recode a categorical or binary variable so that reference level is the mode in the meta-analysis
code_mode_as_ref = function(vec, isMeta){
  
  # # TEST ONLY
  # vec = d$own_mother
  # isMeta = d$isMeta
  
  # make sure the variable is character
  if ( mode(vec) == "logical" ) {
    vec[ vec == TRUE ] = "yes"
    vec[ vec == FALSE ] = "no"
  }
  if( mode(vec) != "character" ) stop( "vec needs to be character or logical" )
  
  # order the character vector 
  correctOrder = names( rev( sort( table( vec[ isMeta == TRUE ] ) ) ) )
  
  # prepend letters to the ordered levels to retain the ordering
  # exploits R's natural alphabetical ordering of factors 
  alphaLevels = paste( letters[ 1:length(correctOrder) ], correctOrder, sep = "." )
  
  # recode values of vec with the letter-prepended versions
  vec2 = vec
  for ( i in 1:length(correctOrder) ) {
    vec2[ vec2 == correctOrder[i] ] = alphaLevels[i]
  }
  
  # sanity check
  #table(vec, vec2)
  
  return(vec2)
}

#code_mode_as_ref(vec = "own_mother", d$isMeta)