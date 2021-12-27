#Important info in Realhelp()
#ALL COMANDS AVAILABLE IN Realtools()

Realhelp <- function(){
  print("for qPCR relative quantification (Pffafl & Livak Method)")
  print("Input csv. dataframe 1 (the qPCR results) columns must have, respectively: target gene Ct, control gene Ct, tissue")
  print("Input csv. dataframe 2 (the standard curve) columns must have, respectively: target gene Ct, control gene Ct, amount of RNA/DNA/cDNA.")
  print("When using Livak method, dataframe 2 it's not needed")
  print("ALL COMMANDS AVAILABLE IN Realtools()")}

Realtools <- function(){
  print("For Pfaffl method use qPCR_1(Df1,Df2,Target gene,Control gene,Target tissue, Control tissue)")
  print("For Livak method use qPCR_2(Df1,Target gene,Control gene,Target tissue, Control tissue)")
  print("To save Pfaffl results in csv. use Save1()")
  print("To save Livak results in csv. use Save2()")
  print("For efficiencies use dataframe 2 and Efi1(Df2) & Efi2(Df2) commands")}

qPCR_1 <- function(a,b,x,y,c,d){
  Re <- data.frame(
    (10^(-1/(unname(round((lm(formula = b[[1]]
                              ~ log10(b[[3]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == c]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == c])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == c]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == c])[["conf.int"]][1:2])),

    (10^(-1/(unname(round((lm(formula = b[[2]]
                              ~ log10(b[[3]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == d]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == d])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == d]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == d])[["conf.int"]][1:2])),


    (10^(-1/(unname(round((lm(formula = b[[1]]
                              ~ log10(b[[3]])))$coeff[2],digits = 2)))))^

      (c(unname((t.test(a[,x][a[,3] == c]))[["estimate"]])
         ,
         t.test(a[,x][a[,3] == c])[["conf.int"]][1:2])-
         c(unname((t.test(a[,y][a[,3] == c]))[["estimate"]])
           ,
           t.test(a[,y][a[,3] == c])[["conf.int"]][1:2]))/
      (10^(-1/(unname(round((lm(formula = b[[2]]
                                ~ log10(b[[3]])))$coeff[2],digits = 2)))))^

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
  View(qPCR_Results_Pfaffl)}

qPCR_2 <- function(a,x,y,c,d){
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

Efi1 <- function (b){10^(-1/(unname(round((lm(formula = b[[1]]
                                              ~ log10(b[[3]])))$coeff[2],digits = 2))))}

Efi2 <- function (b){10^(-1/(unname(round((lm(formula = b[[2]]
                                              ~ log10(b[[3]])))$coeff[2],digits = 2))))}

Save1 <- function(){write.csv(qPCR_Results_Pfaffl,"qPCR_Results_Pfaffl.csv",row.names = FALSE)
}

Save2 <- function(){write.csv(qPCR_Results_Livak,"qPCR_Results_Livak.csv",row.names = FALSE)
}
