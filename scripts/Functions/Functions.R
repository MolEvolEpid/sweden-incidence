# #install dependencies 
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("readr")) install.packages("readr")
if(!require("dplyr")) install.packages("dplyr")
if(!require("readxl")) install.packages("readxl")

#Load infection age distribution information 
infection.ages.gamma.all.new.20231024 <- readRDS("infection.ages.gamma.all.new.20231024.rds")

#Load in IDs for individulas who use a typical distribution
Typical_dist_ID<-readRDS("Typical_dist_ID.rds")

#Add additional variable to infection ages with dates scaled by diagnosis data

for (i in 1:dim(infection.ages.gamma.all.new.20231024)[1]){
  infection.ages.gamma.all.new.20231024$infection.age[[i]]$sx<-infection.ages.gamma.all.new.20231024$diagnosis_date[i]-infection.ages.gamma.all.new.20231024$infection.age[[i]]$x
}


#function to approximate the pdf and cdf for analysis
pdfapprox<-function(dist){
  max_age<-max(dist$sx) #find max value of infection age distribution
  func<-approxfun(dist$sx,dist$y,method="linear", yleft = 0, yright = 0,f=1) #approximate distribution
  integrate_const<-integrate(func,lower=min(dist$sx),upper=max(dist$sx), stop.on.error = FALSE, subdivisions = 10000) #integrate over distribution
  dist$Ny<-dist$y/integrate_const$value #Normalise values to have a total value of 1
  dist$dfunc<-approxfun(dist$sx,dist$Ny,method="linear", yleft = 0, yright = 0,f=1) #approximate function again to get true approximation
  
  store=c() #make storage array to form cumulative function
  
  for (k in dist$sx){
    if (k==min(dist$sx)){
      store=append(store,0)
    }
    else{
      tmp<-integrate(dist$dfunc,lower=min(dist$sx),upper=k, stop.on.error = FALSE, subdivisions = 10000) #calculate CDF for probability of diagnosis
      store=append(store,tmp$value)
    }
  }
  dist$pfunc=approxfun(dist$sx,store,method="linear", yleft = 0, yright = 1,f=1) #Save and approximaate function
  dist$cy<-dist$pfunc(dist$sx)
  
  remove=c() #find irregularrities in the monoticity of a CDF and smoth them out
  for (i in 1:length(dist$sx)){
    
    if (i==1){
      if(dist$cy[i]>dist$cy[i+1]){
        remeove=append(remove,i)
      }
    }
    else if(i==length(dist$sx)){
      if(dist$cy[i-1]>dist$cy[i]){
        
        remeove=append(remove,i)
      }
    }
    else{
      if(!(dist$cy[i-1]<=dist$cy[i] & dist$cy[i]<=dist$cy[i+1])){
        remove=append(remove,i)
      }
    }
  }
  return(dist) #return a function with pdf, and CDF
}

#function to provide posterior probabilities of conditions

#' @param dates A vector containing the dates at which an individual was diagnosed with a condition.
#' @param conditions A vector contain the corresponding conditions to match the order in the dates varaible.
#' @param patID A single integer value for the patient ID as listed in the INFCareHIV data sets.
#' @param plots A Boolean expression for whether plots are required, default is FALSE.
#' @param csv A Boolean expression for whether a .csv file is required, default is FALSE. Otherwise a dataframe will be printed
#' @export
#' 
#' 
Condprobs<-function(dates,conditions,patID,plots=FALSE,csv=FALSE){

  #indicator for typical distribution
  typical=FALSE
  
  #Find patient ID in infection age distributions 
  ID<-which(infection.ages.gamma.all.new.20231024$patient_ID==patID) 
  
  #check if using a typical distribution
  if (ID %in% Typical_dist_ID){
    typical=TRUE
  }
  
  #Approximate the pdf function
  pdf<-pdfapprox(infection.ages.gamma.all.new.20231024$infection.age[[ID]]) 
  
  #pull diagnosis date
  diag<-infection.ages.gamma.all.new.20231024$diagnosis_date[ID] 


    df <- data.frame(ID=double(),
                     Date=as.Date(character()),
                     Condition=character(), 
                     Probability=double(), 
                     NO_Probability=double(),
                     Typical=character(),
                     stringsAsFactors=FALSE) #empty storage array to show condition probabilities with typical
  
  
  #convert condition dates into correct format
  cdates<-as.double(as.Date(dates, "%Y-%m-%d"))/365.25+1970 
  
  #loop over condition dates
  for(i in 1:length(cdates)){
    probs<-round(integrate(pdf$dfunc,lower = min(pdf$sx), upper = cdates[i], stop.on.error = FALSE, subdivisions = 10000)$value,3)
    if(typical==FALSE){
      df<-rbind(df,c(ID,dates[i],conditions[i],probs,1-probs,"FALSE"))
    }
    else{
      df<-rbind(df,c(ID,dates[i],conditions[i],probs,1-probs,"TRUE"))
    }
  }
  


    colnames(df)[1]<-"ID"
    colnames(df)[2]<-"Dates"
    colnames(df)[3]<-"Condition"
    colnames(df)[4]<-"Posterior probability of HIV infection before diagnosis"
    colnames(df)[5]<-"Posterior probability of no HIV infection before diagnosis"
    colnames(df)[6]<-"Typical used"

    if(plots==TRUE){
      plotdf<-data.frame(x=pdf$sx,y=pdf$dfunc(pdf$sx))
      
      pdf( paste0("Comorbidities_pat_",ID,".pdf"), width = 11, height = 8 )
      
      for (i in 1:length(cdates)){
        
        #plot figures with individual probabilities
        p<-ggplot(plotdf, aes(x = x, y = y)) + 
          geom_line(stat = 'identity') +
          geom_ribbon(aes(ymin=0,ymax=if_else(x<=cdates[i],y,NA),fill="#0BB702"))+
          geom_ribbon(aes(ymin=0,ymax=if_else(x>cdates[i],y,NA),fill="#F8766D"))+
          geom_vline(xintercept = cdates[i])+ 
          guides(fill=guide_legend(title="Probability"))+ 
          xlab("Year") +
          ylab("Density")+
          ggtitle(paste("Probability of HIV infection at diagnosis of",df[i,3],"on the",df[i,2]))+ 
          theme_classic()+
          scale_fill_manual(values=c("#0BB702", "#F8766D"), labels=c(paste(df[i,4],"Yes"),paste(df[i,5],"No")),name="Probs")
        
        print(p)
      }
      
      dev.off()
    }
  
  if(csv==TRUE){
    return(df)
  }
  else{
    print(df)
  }
  
  

  
}  

all_pat<-function(Pat_conditions,plots=FALSE,csv=TRUE){
  Pat_conditions$Date<-as.character(Pat_conditions$Date)
  #find unique patient IDS
  Pat_IDS<-Pat_conditions %>% distinct(INFCare_ID)
  for (j in 1:dim(Pat_IDS)[1]){
    temp<-Pat_conditions %>% filter(INFCare_ID==Pat_IDS$INFCare_ID[j])
    out<-suppressWarnings(Condprobs(temp$Date,temp$Condition,Pat_IDS$INFCare_ID[j],plots=plots,csv=csv))
    if(j==1){
      df<-out
    }
    else df<-rbind(df,out)
  }
  write_csv(df,"Comorbidities_out.csv")
}

Single_pat<-function(ID,Pat_conditions,plots=FALSE,csv=FALSE){
  temp<-Pat_conditions %>% filter(INFCare_ID==ID)
  Condprobs(as.character(temp$Date),temp$Condition,ID,plots,csv)
}






