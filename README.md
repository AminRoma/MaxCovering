---
title: "Maximal Coverage Problem Optimization Code"
output: html_document
author: "Amin Asadi"
---  

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Packages

This is a document to illustrate our code to solve **maximal coverage for a capacitated facility location problem**. First we load the packages and note that for the first time use, we need to install packages using `install.packages()` command. 


```{r packages, eval = FALSE}
library(listviewer)
library (Rglpk)
library(dplyr)
library(magrittr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)
library(gtools)
library(R.utils)
library(reader)
library(stringr)
```
## Paramters

We have 42 candidate facility location with three level of services (1- buliding no service, 2- partial service, 3- Full service). Besides, we have 50 existing facilities that are extendable. Hence, we have 176 candidate in total. The number of exisiting facilities are 168. We assume the overload of a parking facility can be covered by another facility within the `max_Distance` miles. The available budget for the expansion is `max_budget` dollars.  


```{r parameters, eval=FALSE}
rm(list=ls())
M1 <-176
N1 <-168
max_Distance<-40
max_budget<-10000000
Start = Sys.time()
```
### Reading from the simulation or adjusting based on the simulation
```{r Simulation, eval=FALSE}
Tick<-10
Number_Tick_Per_Hours <-60/Tick
Simulation_Day<-14
Simulation_Hours<-Simulation_Day*24
Simulation_Lenght<-Simulation_Hours*Number_Tick_Per_Hours
# if report name Report=0 and if we use IDs report=1.
Report<-1
```
### Read distance matrix, cost and capacity vectors and station results

```{r data, eval=FALSE}
setwd("/Users/asadi/Desktop/MaxCoveringallFiles")
Sys.setenv('R_MAX_VSIZE'=64000000000)
Sys.getenv('R_MAX_VSIZE')
Distances <- read.delim("DistanceFull.txt", sep = "", header= FALSE)
CostsCapacities <- read.delim("CostCapacity.txt", sep = "", header= FALSE)
D = read.csv("stationsresultsSep4.csv", header =TRUE)
Facilities_176 <- read.delim("176Facilities.txt", sep = "\t", header= FALSE)
Facilities_168 <- read.delim("168Facilities.txt", sep = "\t", header= FALSE)
```

### Create build and expnasion cost `f` and capacity of candidates `q`

```{r vectors, eval=FALSE}
f <-CostsCapacities$V1
f <-as.matrix(f)
q <-CostsCapacities$V2
q <-as.matrix(q)
```

### Make the correct order of facilities and Create overcapacity Matrix
```{r dj, eval=FALSE}
D = D[mixedorder(as.character(D$Parking.ID)),]
D$ID <- seq.int(nrow(D))
D$ID <-(D$ID %/%Simulation_Lenght)+1
D$OverCapacity<- with(D, round(pmax(Ratio-1,0)*Parking.Capacity))

h<-runif(N1, 0)
d<-runif(N1, 0)

# d[k] is the vector of overcapcity facilities
for (k in 1:N1){
  h[k]<-D %>%
    filter(ID==k) %>%
    select(OverCapacity) %>%
    summarise(x = max(OverCapacity))
  d[k]<-h[[k]]
}

```
### Help solver
We remove the exisiting facilities with zero demand (overcapacity) from `d` and associted distances from `Distance` matrix and save them into `d1` and `Distance1`, respectively. We keep the index of remaining nodes into `NamesVec`, and use it when are reporting our solution.

```{r Helpsolver, eval=FALSE}
d1 <-data.frame()[1,]
Distances1<-data.frame()[1:M1, ]
NamesVec <-data.frame()[1,]
for (j in 1:N1){
  if (d[j]!=0){
    Distances1<-cbind(Distances1, Distances[,j])
    d1<-cbind(d1,d[j])
    NamesVec<-cbind(NamesVec,j)
  }
}

NamesVec<-t(NamesVec)
d1<-t(d1)
rownames(d1)<-c(NamesVec)

rownames(Distances1)<-(1:M1)
colnames(Distances1)<-c(NamesVec)

# N2 is the number of exisiting facilities with demand not equal to 0.
N2<-as.numeric(ncol(Distances1))

# w(i,j) is the distance between nodes i and j in the MIP model.
w<-function(i,j){
  vapply(seq_along(i), function(k) Distances1[i[k], j[k]], numeric(1L))
}

```
### MIP model

We create our model using `ompr` package and solve it with `glpk` solver. We use `sink()` command to create a log file for the solver report.

```{r MIP, eval=FALSE}
sink("./logofcode.txt")
model<-MIPModel() %>%
  add_variable(y[i], i = 1:M1, type = "binary") %>%
  add_variable(x[i,j], i = 1:M1, j=1:N2, type = "continuous", lb = 0, ub =1 ) %>%
  add_variable(z[i,j], i = 1:M1, j=1:N2, type = "binary") %>%

  set_objective(sum_expr(x[i,j]*d1[j], i = 1:M1, j=1:N2)/sum_expr(d1[j], j=1:N2), "max") %>%

  add_constraint(z[i,j]*w(i,j) <= y[i] * max_Distance , j=1:N2, i = 1:M1) %>%
  add_constraint(x[i,j]-z[i,j] <=0.0  , i = 1:M1, j=1:N2) %>%
  add_constraint(sum_expr(f[i]*y[i], i = 1:M1)<= max_budget) %>%
  add_constraint(sum_expr(x[i, j], i = 1:M1) <= 1.0, j = 1:N2) %>%
  add_constraint(sum_expr(d1[j]*x[i, j], j = 1:N2)<= y[i] * q[i] , i = 1:M1) %>%
  add_constraint(y[k]+y[k+42]+y[k+84] <= 1.0, k = 1:42) 
model
# Time limit is 3600 second or an hour.
result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE, tm_limit = 1000*100))
sink()
```
### Refine solver report for near optimal solution

We report the optimality gap for the near optimal solution using the string saved in `OptimalityGap3`.

```{r logfile, eval=FALSE}
LogfileLines<-as.numeric(countLines("logofcode.txt", chunkSize=5e+07))
my_txt_ex2 <- n.readLines(paste("logofcode.txt", sep = ""),header = FALSE, n = 1,  skip =LogfileLines-3 )
startPercentage<-str_locate(my_txt_ex2,"%")[1,1]
endPercentage<-str_locate(my_txt_ex2,"%")[1,2]
OptimalityGap<-substring(my_txt_ex2, startPercentage-4, endPercentage)
OptimalityGap3<-paste("near optimal with optimality gap of", OptimalityGap)

if (result$objective_value >0 & result$objective_value <1 & result[["status"]] =="infeasible"){
  result[["status"]] = OptimalityGap3
}else if (result[["status"]] =="infeasible" & max_budget>1000000){
  result[["status"]] = "opt"
  cat("reduce maximum budget from $", max_budget, "\n")
}
````

### Debug solver for very high maximum budget and unused resources

We observe that when the budget is very high and an amount of budget remain unused while the demand is satisfied with a significant lower budget, the solver report infeasible. In this section we reduce budget behind scence to find a solution and avoid ompr or glpk bug.

```{r logfile2, eval=FALSE}

while(result[["status"]]=="opt"){
  
  rm(model)
  rm(result)
  
  max_budget = max_budget - 2000000

  sink("./logofcode.txt")
  model<-MIPModel() %>%
    add_variable(y[i], i = 1:M1, type = "binary") %>%
    add_variable(x[i,j], i = 1:M1, j=1:N2, type = "continuous", lb = 0, ub =1 ) %>%
    add_variable(z[i,j], i = 1:M1, j=1:N2, type = "binary") %>%
    
    set_objective(sum_expr(x[i,j]*d1[j], i = 1:M1, j=1:N2)/sum_expr(d1[j], j=1:N2), "max") %>%
    
    add_constraint(z[i,j]*w(i,j) <= y[i] * max_Distance , j=1:N2, i = 1:M1) %>%
    add_constraint(x[i,j]-z[i,j] <= 0  , i = 1:M1, j=1:N2) %>%
    add_constraint(sum_expr(f[i]*y[i], i = 1:M1)<= max_budget) %>%
    add_constraint(sum_expr(x[i, j], i = 1:M1) <= 1 , j = 1:N2) %>%
    add_constraint(sum_expr(d1[j]*x[i, j], j = 1:N2)<= y[i] * q[i] , i = 1:M1) %>%
    add_constraint(y[k]+y[k+42]+y[k+84] <= 1.0, k = 1:42) 
  model
  result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE, tm_limit = 1000*300))
  sink()
  
  # Create log file
  LogfileLines<-as.numeric(countLines("logofcode.txt", chunkSize=5e+07))
  my_txt_ex2 <- n.readLines(paste("logofcode.txt", sep = ""),header = FALSE, n = 1,  skip =LogfileLines-3 )
  startPercentage<-str_locate(my_txt_ex2,"%")[1,1]
  endPercentage<-str_locate(my_txt_ex2,"%")[1,2]
  OptimalityGap<-substring(my_txt_ex2, startPercentage-4, endPercentage)
  OptimalityGap3<-paste("near optimal with optimality gap of", OptimalityGap)
  
  if (result$objective_value >0 & result$objective_value <1 & result[["status"]] =="infeasible"){
    result[["status"]] = OptimalityGap3
  }else if (result[["status"]] =="infeasible" & max_budget>1000000){
    result[["status"]] = "opt"
    cat("reduce maximum budget from $", max_budget, "\n")
  }
}

```

###  Make data ready for presenting solution

```{r solution, eval=FALSE}

matching <- result 
soln.x <- get_solution(result, x[i, j])
soln.y <- get_solution(result, y[i])
lapply(matching$solution,head)
summary(matching)

# Debug the solver report that shows infeasible when it can not find the optimal while it is feasible
result
  if (result$objective_value >0 & result$objective_value <1 & result[["status"]] =="infeasible")
  result[["status"]] = OptimalityGap3

# write x[i,j] and y[i] solutions to .txt files
soln.x <- get_solution(result, x[i, j])
write.table(soln.x, "Xij.txt", sep="\t")
soln.y <- get_solution(result, y[i])
write.table(soln.y, "Yi.txt", sep="\t")



# Find facilities name from 168 facilities based on what we have in the demand Vec
DemandFacWithNames<-Facilities_168[Facilities_168$V1 %in% NamesVec,]
DemandFacWithNames<-cbind(DemandFacWithNames, Temp_ID = seq(1:nrow(NamesVec)))

# Map solutions with ID and names
soln.y<-cbind(soln.y, Facilities_176$V2)
soln.y<-rename(soln.y, "ID" = "Facilities_176$V2")
soln.y$ID<-as.character(soln.y$ID)
soln.y<-cbind(soln.y, Facilities_176$V3)
soln.y<-rename(soln.y, "Name" = "Facilities_176$V3")
soln.y$Name<-as.character(soln.y$Name)

# Divide candiates based on 4 level of service (y1:No serive, y2: Partial Service, y3: Full service, y4: Expansion)
y1<-filter(soln.y, soln.y$i <  43  & soln.y$value==1)
y2<-filter(soln.y, soln.y$i >= 43  & soln.y$i<=85 & soln.y$value==1)
y3<-filter(soln.y, soln.y$i >= 85  & soln.y$i<=126 & soln.y$value==1)
y4<-filter(soln.y, soln.y$i >= 127 & soln.y$value==1)

# Convert x[i,j] to x[name/ID of supply, name/ID of demand]
ID_Supplies<-as.character(rep(Facilities_176$V2, times = N2))
ID_Demands<-as.character(rep(DemandFacWithNames$V2, each = M1))
Name_Supplies<-as.character(rep(Facilities_176$V3, times = N2))
Name_Demands<-as.character(rep(DemandFacWithNames$V3, each = M1))
soln.x<-cbind(soln.x, ID_Supplies,ID_Demands, Name_Supplies, Name_Demands)


x1<-filter(soln.x, soln.x$i <  43  & soln.x$value!=0)
x2<-filter(soln.x, soln.x$i >= 43  & soln.x$i<=85 & soln.x$value!=0)
x3<-filter(soln.x, soln.x$i >= 85  & soln.x$i<=126 & soln.x$value!=0)
x4<-filter(soln.x, soln.x$i >= 127 & soln.x$value!=0)
```

###  Create report

```{r report, eval=FALSE}

# 0 if report based on name and 1 if we report based on ID
Report<-1
myFunction<-function(result, Report){
  
  if (result[["status"]] == "infeasible"){
    cat("No solution is found with the current budget. Please increase the budget to find a feasible solution.", "\n")
  } else if(result[["status"]] == "opt"){
    cat("This problem has optimal solution. Please reduce the budget to see the optimal solution.", "\n")
  }else{
    
    if (nrow(y1)==0 & nrow(y2)==0 & nrow(y3)==0 & nrow(y4)==0 )
      print("No change is necessary")
    else{
      
      if (Report ==1){
        if (nrow(y1)>0)
          for (j in 1:nrow(y1))
            cat("Build facility", y1[j,4], "with no service", "\n")
        if (nrow(y2)>0)
          for (j in 1:nrow(y2))
            cat("Build facility", y2[j,4], "with partial service", "\n")
        if (nrow(y3)>0)
          for (j in 1:nrow(y3))
            cat("Build facility", y3[j,4], "with full service", "\n")
        if (nrow(y4)>0)
          for (j in 1:nrow(y4))
            cat("Expand facility", y4[j,4], "\n")
        
        
        
        if (nrow(x1)>0)
          for (j in 1:nrow(x1))
            cat(round(x1[j,4]*100,digits = 2),"% of demand for facility", as.character(x1[j,6]), "is satisfied by buliding new facility",  as.character(x1[j,5])," with no serivce","\n")
        if (nrow(x2)>0)
          for (j in 1:nrow(x2))
            cat(round(x2[j,4]*100,digits = 2),"% of demand for facility", as.character(x2[j,6]), "is satisfied by buliding new facility",  as.character(x2[j,5])," with partial serivce", "\n")
        if (nrow(x3)>0)
          for (j in 1:nrow(x3))
            cat(round(x3[j,4]*100,digits = 2),"% of demand for facility", as.character(x3[j,6]), "is satisfied by buliding new facility",  as.character(x3[j,5])," with full serivce","\n")
        if (nrow(x4)>0)
          for (j in 1:nrow(x4))
            cat(round(x4[j,4]*100,digits = 2),"% of demand for facility", as.character(x4[j,6]), "is satisfied by expanding facility",  as.character(x4[j,5]), "\n")
      }
      
      
      else {
        if (nrow(y1)>0)
          for (j in 1:nrow(y1))
            cat("Build facility", y1[j,5], "with no service", "\n")
        if (nrow(y2)>0)
          for (j in 1:nrow(y2))
            cat("Build facility", y2[j,5], "with partial service", "\n")
        if (nrow(y3)>0)
          for (j in 1:nrow(y3))
            cat("Build facility", y3[j,5], "with full service", "\n")
        if (nrow(y4)>0)
          for (j in 1:nrow(y4))
            cat("Expand facility", y4[j,5], "\n")
        
        
        
        if (nrow(x1)>0)
          for (j in 1:nrow(x1))
            cat(round(x1[j,4]*100,digits = 2),"% of demand for facility", as.character(x1[j,8]), "is satisfied by buliding new facility",  as.character(x1[j,7])," with no serivce","\n")
        if (nrow(x2)>0)
          for (j in 1:nrow(x2))
            cat(round(x2[j,4]*100,digits = 2),"% of demand for facility", as.character(x2[j,8]), "is satisfied by buliding new facility",  as.character(x2[j,7])," with partial serivce", "\n")
        if (nrow(x3)>0)
          for (j in 1:nrow(x3))
            cat(round(x3[j,4]*100,digits = 2),"% of demand for facility", as.character(x3[j,8]), "is satisfied by buliding new facility",  as.character(x3[j,7])," with full serivce","\n")
        if (nrow(x4)>0)
          for (j in 1:nrow(x4))
            cat(round(x4[j,4]*100,digits = 2),"% of demand for facility", as.character(x4[j,8]), "is satisfied by expanding facility",  as.character(x4[j,7]), "\n")
        
      }
    }
  }
  
  a = 0;
  
  for (i in 1:M1) {
    a = a + f[i]*soln.y[i,3]
  }
  
  if (result[["status"]] != "opt"){
    cat("Budget to spend on the new and expanded facilities is $",a, "\n")
    cat("This solution is", result$status, "\n")
    cat("The average met demand is", round(result$objective_value*100,digits = 2),"%", "\n")
  } else   
  {
    cat("The average percentage of met demand with this budget is 100%", "\n")
  }
}


myFunction(result, Report)


End = Sys.time()
End - Start 
```
