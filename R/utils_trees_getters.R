# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions to query the properties of a tree, to check its
# consistency etc. All these functions assume to be working
# with the adjacency matrix representation of a tree
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Return true if M is a tree
is_tree = function(M, allow.empty = FALSE)
{
  if(sum(M) == 0 & allow.empty) return(TRUE)
  
  cS = colSums(M)
  
  # - No more than 1 parent
  # - 1 root
  # - No nodes with 0 parents
  (max(cS) == 1) && (length(root(M)) == 1) && (sum(cS) == length(cS) - 1)
}

# Compute the root(s) of a model
root = function(model){
  s = colSums(model)
  return(names(s[s==0]))
}

# Compute the children of "var" in a model
children = function(model, var)
{
  model = model[var, ]
  model = model[model == 1]
  
  if(is.null(model)) return(NULL)
  return(names(model))
}

# Return the parent of "variable" in model
pi = function(model, variable)
{
  model = model[, variable]
  if(any(model > 0)) model = model[model > 0]
  else return(NULL)
  return(names(model))
}

# Compute the leaves of a model
leaves = function(model){
  s = rowSums(model)
  model = t(model)
  return(names(s[s==0]))
}

# Given a Data Frame representation of a model, compute the set of
# nodes reachable from x
reach = function(df, x)
{
  if(!any(df$from == x)) return(NULL)
  
  dfB = df[ df$from == x, , drop = F]
  r = dfB$to
  
  # print(dfB)
  for(i in 1:length(r))
  {
    r = c(r, reach(df, r[i]))
  }
  
  return(r)
}
