Realhelp <- function(){
  print("for qPCR relative quantification (Pffafl & Livak Method)")
  print("Input csv. dataframe 1 (the qPCR results) columns must have, respectively: target gene Ct, control gene Ct, tissue")
  print("Input csv. dataframe 2 (the standard curve) columns must have, respectively: target gene Ct, control gene Ct, amount of RNA/DNA/cDNA.")
  print("When using Livak method, dataframe 2 it's not needed")}

qPCR_1 <- function(a,b,x,y,c,d){
  install.packages("dplyr")
  library("dplyr")
  Re <- data.frame(
    (10^(-1/(unname(round((lm(formula = as.list(b %>% select(1))[[1]]
                              ~ log10(as.list(b %>% select(3))[[1]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == c]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == c])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == c]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == c])[["conf.int"]][1:2])),

    (10^(-1/(unname(round((lm(formula = as.list(b %>% select(2))[[1]]
                              ~ log10(as.list(b %>% select(3))[[1]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == d]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == d])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == d]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == d])[["conf.int"]][1:2])),


    (10^(-1/(unname(round((lm(formula = as.list(b %>% select(1))[[1]]
                              ~ log10(as.list(b %>% select(3))[[1]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == c]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == c])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == c]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == c])[["conf.int"]][1:2]))/
      (10^(-1/(unname(round((lm(formula = as.list(b %>% select(2))[[1]]
                                ~ log10(as.list(b %>% select(3))[[1]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == d]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == d])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == d]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == d])[["conf.int"]][1:2]))

  )
  colnames(Re)<-c("Target_tissue","Control_tissue","Ratio")
  rownames(Re)<-c("Mean","Lower CI","Upper CI")
  qPCR_Results_Pfaffl <<- Re
  View(qPCR_Results_Pfaffl)
}

qPCR_2 <- function(a,x,y,c,d){
  install.packages("dplyr")
  library("dplyr")
  Re <- data.frame(
    (2)^

      (c(unname((t.test(a[,x][a[,3] == c]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == c])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == c]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == c])[["conf.int"]][1:2])),

    (2)^

      (c(unname((t.test(a[,x][a[,3] == d]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == d])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == d]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == d])[["conf.int"]][1:2])),


    (2)^

      (c(unname((t.test(a[,x][a[,3] == c]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == c])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == c]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == c])[["conf.int"]][1:2]))/
      (2)^

      (c(unname((t.test(a[,x][a[,3] == d]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == d])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == d]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == d])[["conf.int"]][1:2]))

  )
  colnames(Re)<-c("Target_tissue","Control_tissue","Ratio")
  rownames(Re)<-c("Mean","Lower CI","Upper CI")
  qPCR_Results_Livak <<- Re
  View(qPCR_Results_Livak)
}

Save <- function(){
  write.csv(qPCR_Results_Pfaffl,"qPCR_Results_Pfaffl.csv")
  write.csv(qPCR_Results_Livak,"qPCR_Results_Livak.csv")
}
