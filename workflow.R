
path_input <- "a_runFolder"
path_input_occ <- "a_runFolder/input_occ"
path_input_index <- "a_runFolder/index"

dir.create(path_input)
dir.create(path_input_occ)
dir.create(path_input_index)
GLOBAL_variables =  paste0("bio",c(10,11,16,17) ) 
GLOBAL_resolution = 2.5
GLOBAL_bgnum = 10000
wc_raw <- paste0("data/2_wc/wc1/",GLOBAL_resolution,"/",GLOBAL_variables,".bil")
wc_raw <- stack(wc_raw)

all_spp <- list.files("data/3_soberon2017/occurrences/locs_rem/")
all_sppcsv <- list.files("data/3_soberon2017/occurrences/locs_rem/",full.names = T)
length(all_sppcsv)
has_M <- (file.exists(paste0("data/3_soberon2017/M/",gsub(" .csv",".shp",all_spp))))
table(has_M) 
all_spp <- all_spp[has_M]
all_sppcsv <- all_sppcsv[has_M]
all_spp_sp <- paste0(path_input_occ,"/",basename(all_sppcsv),".rds" )
all_M <- paste0("data/3_soberon2017/M/",gsub(" .csv",".shp",all_spp))
occ_summary <- data.frame(species= gsub(" .csv","",all_spp),
                          csv_occ= all_sppcsv,
                          occ_sp= all_spp_sp,
                          path_M = all_M,
                          num_occ= 0,
                          num_occ_clean = 0,
                          stringsAsFactors = F)
sss=1
for (sss in 1: nrow(occ_summary) ){
  print(sss)
  myocc <- read.csv(occ_summary$csv_occ[sss])  
  coordinates(myocc) <-~ long + lat  
  myocc_sp <- extract(wc_raw,myocc,cellnumbers=TRUE,sp=TRUE)
  dups <- duplicated(myocc_sp$cells)
  myocc_sp <- myocc_sp[!dups,]
  myocc_sp <- myocc_sp[!is.na(myocc_sp$bio10),]
  occ_summary$num_occ[sss] <- nrow(myocc)
  occ_summary$num_occ_clean[sss] <- nrow(myocc_sp)  
  saveRDS(myocc_sp,occ_summary$occ_sp[sss])
}
saveRDS(occ_summary,paste0(path_input_index,"/occ_summary_v0.rds"))
nrow(occ_summary) 
occ_summary <- subset(occ_summary,num_occ_clean>=100)
bad_species <- c("9323_Anthracothorax_dominicus",
                 "9426_Chlorostilbon_ricordii",
                 "9427_Chlorostilbon_swainsonii",
                 "10056_Calliphlox_evelynae",
                 "10067_Mellisuga_minima",
                 "9327_Anthracothorax_mango",
                 "9340_Eulampis_holosericeus",
                 "9346_Orthorhynchus_cristatus"
)
occ_summary <- subset(occ_summary, ! species%in% bad_species )
saveRDS(occ_summary,paste0(path_input_index,"/occ_summary_v1.rds"))
nrow(occ_summary)

myWorkflow <- function(this_species, 
                       flag_clamp, 
                       flag_testP="testOnly",
                       flag_testA,  
                       path_input, 
                       buffersize_all =c(-9999,
                                         10,15,20,30,500,800,1500,
                                         50,100,200
                       ) ,
                       cut_all=c("block",
                                 "random"
                       ), 
                       rep_all = 4,
                       MYFEATURES = "lqh",
                       buffer_method = "trainingOcc", 
                       all_algorithm= c("maxnet"), 
                       cup_num = 7,
                       flag_redo_level = 10 
){  

catf <- function(..., file="temp/temp_preparebuffer_v1.log", append=TRUE){
  cat(..., file=file, append=append)
}
my_blockcut<- function (occ) {
  occ <- as.data.frame(occ)
  rownames(occ) <- 1:nrow(occ)
  noccs <- nrow(occ)
  n1 <- ceiling(nrow(occ)/2)
  n2 <- floor(nrow(occ)/2)
  n3 <- ceiling(n1/2)
  n4 <- ceiling(n2/2)
  grpA <- occ[order(occ[, 2]), ][1:n1, ]
  grpB <- occ[rev(order(occ[, 2])), ][1:n2, ]
  grp1 <- grpA[order(grpA[, 1]), ][1:(n3), ]
  grp2 <- grpA[!rownames(grpA) %in% rownames(grp1), ]
  grp3 <- grpB[order(grpB[, 1]), ][1:(n4), ]
  grp4 <- grpB[!rownames(grpB) %in% rownames(grp3), ]
  bvert <- mean(c(max(grp1[, 1]), min(grp2[, 1])))
  tvert <- mean(c(max(grp3[, 1]), min(grp4[, 1])))
  horz <- mean(c(max(grpA[, 2]), min(grpB[, 2])))
  
  r <- data.frame()
  if (nrow(grp1) > 0) 
    grp1$grp <- 1
  r <- rbind(r, grp1)
  if (nrow(grp2) > 0) 
    grp2$grp <- 2
  r <- rbind(r, grp2)
  if (nrow(grp3) > 0) 
    grp3$grp <- 3
  r <- rbind(r, grp3)
  if (nrow(grp4) > 0) 
    grp4$grp <- 4
  r <- rbind(r, grp4)
  occ.grp <- r[order(as.numeric(rownames(r))), ]$grp
  out <- list(occ.grp = occ.grp)
  return(out)
}
calPAUC <- function (p,a,model,E=0.95,flag_clamp=FALSE,type){
  p <- predict( object= model,newdata=p,type=type, clamp=flag_clamp)
  a <- predict( object= model,newdata=a,type=type, clamp=flag_clamp)
  input1 <- c ( rep(1,length(p) ),rep(0,length(a) ) )
  input2 <- c(p,a)    
  auc_test_p = tryCatch({ pROC::auc(input1, input2,
                                    partial.auc=c(1, E), partial.auc.focus="se",
                                    partial.auc.correct=T)[1]
  }, error=function(e)e )  
  if(! is.numeric(auc_test_p) ) {
    auc_test_p= NA
    }  
  return(auc_test_p)
}
get.params <- function(model){
  lambdas <- model@lambdas[1:(length(model@lambdas)-4)]
  countNonZeroParams <- function(x) {
    if (strsplit(x, split=", ")[[1]][2] != '0.0') 1
  }
  no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(no.params)
}
calc.aicc <- function(model, statBG,nocc,training_presences,flag_clamp=TRUE,type){
  nparam <- length(model$betas)  
  predictive.maps <- predict(object=model,
                             newdata=statBG,
                             type=type,
                             clamp=flag_clamp  )  
  probsum <- sum(predictive.maps)  
  vals <- predict(object=model,
                  newdata=training_presences,
                  type=type,
                  clamp=flag_clamp)
  AIC.valid <- nparam < nocc
  LL <- colSums(log(t(t(vals)/probsum)), na.rm=T)
  AICc <- (2*nparam - 2*LL) + (2*(nparam)*(nparam+1)/(nocc-nparam-1))
  AICc[AIC.valid==FALSE] <- NA
  AICc[is.infinite(AICc)] <- NA
  return (AICc)
}
myENM_auc <- function(one_model,occ_train, occ_test,
                      bg_train,
                      bg_test,
                      athreshold=NA,flag_clamp,
                      which_algorothm
                      ){
  if(which_algorothm=="glm"){  
    ped_type="link"
  }
  if(which_algorothm=="maxnet"){
    ped_type="exponent"
  }  
  ev_train <- dismo::evaluate(p =occ_train,a =
                                bg_train,model = one_model, 
                              clamp=flag_clamp  )
  pp <- predict(object=one_model,
                newdata=occ_train,
                clamp=flag_clamp)
  pp <- na.omit(pp)
  ev_MTP <- min(pp)
  ev_5 <- quantile(pp,5/100)
  ev_test_MTP <- dismo::evaluate(p =occ_test,a =bg_test,model = one_model ,tr=ev_MTP,
                                 clamp=flag_clamp)
  one_sen <- ev_test_MTP@TPR
  one_spec <- ev_test_MTP@TNR
  one_TSS <- one_sen + one_spec - 1
  one_kappa <- ev_test_MTP@kappa
  one_AUCtest <- ev_test_MTP@auc
  one_AUCtrain <- ev_train@auc
  one_AUCdiff <- one_AUCtrain - one_AUCtest
  ev_test_5 <- dismo::evaluate(p =occ_test,a =bg_test,model = one_model ,tr=ev_5,
                               clamp=flag_clamp)
  one_TPR_5 <- ev_test_5@TPR
  one_FPR_5 <- ev_test_5@FPR
  one_TNR_5 <- ev_test_5@TNR
  one_FNR_5 <- ev_test_5@FNR
  one_soren_5 <- 2*one_TPR_5 / (2*one_TPR_5 + one_FPR_5 + one_FNR_5)
  
  one_sen_5 <- ev_test_5@TPR
  one_spec_5 <- ev_test_5@TNR
  one_TSS_5 <- one_sen_5 + one_spec_5 - 1
  one_kappa_5 <- ev_test_5@kappa
  pAUC5_train <- calPAUC (p =occ_train,a =bg_train, 
                          model = one_model, E=0.95,
                          flag_clamp=flag_clamp,
                          type=ped_type)
  pAUC5_test  <- calPAUC (p =occ_test, a =bg_test,
                          model = one_model, E=0.95,
                          flag_clamp=flag_clamp,
                          type=ped_type)
  pAUC5_diff <- pAUC5_train-pAUC5_test
  ped_occ_train <- predict(object=one_model,
                          newdata=occ_train,
                          clamp=flag_clamp)
  ped_occ_test <- predict(object=one_model,
                          newdata=occ_test,
                          clamp=flag_clamp)
  ped_bg_test <- predict(object=one_model,
                         newdata=bg_test,
                         clamp=flag_clamp)
  boyce_train = ecospat::ecospat.boyce(fit=c(ped_bg_test,ped_occ_train),obs=ped_occ_train)
  boyce_test = ecospat::ecospat.boyce(fit=c(ped_bg_test,ped_occ_test),obs=ped_occ_test)  
  if(which_algorothm=="maxnet"){
    # AIC
    #npp <- get.params(one_model)
    npp <- length(one_model$betas)
    aic <- calc.aicc(one_model,bg_train,nrow(occ_train),
                     # this index was calculated wrong previously. testing bg was used in previous runs. Training bg should be used here.
                     training_presences=occ_train,
                     flag_clamp=flag_clamp,
                     type=ped_type)
  } else{
    npp = NA
    aic = AIC(one_model)
  }  
  bg_sum_raw <- sum( predict(one_model,bg_train,
                             type=ped_type,
                             clamp=flag_clamp) )  
  one_line <- cbind(
    "n_train_occ"=nrow(occ_train),
    "n_test_occ"=nrow(occ_test),    
    "n_bg"=nrow(bg_train),
    "bg_sum_raw"=bg_sum_raw,
    "sensitivity"=one_sen,
    "specificity"=one_spec,
    "TSS"=one_TSS,
    "Kappa"=one_kappa,
    "sensitivity5"=one_sen_5,
    "specificity5"=one_spec_5,
    "TSS5"=one_TSS_5,
    "Kappa5"=one_kappa_5,
    "AUCtest"=one_AUCtest,
    "AUCtrain"=one_AUCtrain,
    "AUCdiff"=one_AUCdiff,
    "pAUCtest"=pAUC5_test,
    "pAUCtrain"=pAUC5_train,
    "pAUCdiff"=pAUC5_diff,
    "npap"=npp,
    "aic"=aic,
    "one_TPR_5" =one_TPR_5,
    "one_FPR_5" =one_FPR_5,
    "one_TNR_5" =one_TNR_5,
    "one_FNR_5" =one_FNR_5,
    "one_soren_5" =one_soren_5,
    "boyce_test"=boyce_test$cor,
    "boyce_train"=boyce_train$cor
  )
  return(one_line)
}
    wc_raw <- readRDS("data/gis/wc_raw_noisland.rds")   
    GLOBAL_bgnum = 10000
    path_input_occ <- paste0(path_input,"/input_occ")
    path_input_occ_cut <- paste0(path_input,"/input_occ_cut")
    path_input_occ_buff <- paste0(path_input,"/input_occ_buff")
    path_input_map <- paste0(path_input,"/input_map")
    path_input_index <- paste0(path_input,"/index")
    path_output_model <- paste0(path_input,"/output_model")
    path_output_eval <- paste0(path_input,"/output_eval")    
    dir.create(path_input)
    dir.create(path_input_occ)
    dir.create(path_input_index)
    dir.create(path_input_occ_cut)
    dir.create(path_input_occ_buff)
    dir.create(path_input_map)
    dir.create(path_output_model)
    dir.create(path_output_eval)    
    path_temp_raster <- "D:/fx_temp/"
    rasterOptions(tmpdir = path_temp_raster)
    MAXBUFFER <- max(buffersize_all)
  cut_label="block"
  for(cut_label in cut_all){
    cat(cut_label)
    
    myocc_sp <- readRDS(this_species$occ_sp)
    crs(myocc_sp) = crs(wc_raw)
    temp_check_occ = extract(wc_raw[[1]],myocc_sp)
    to_rm_occ = which(is.na(temp_check_occ))
    if(length( to_rm_occ)>0 ){
      print("remove some occ because of NA env value")
      myocc_sp = myocc_sp[-to_rm_occ,]
    }
    path_temp_occ_cut <- paste0(path_input_occ_cut,"/",
                                basename(this_species$csv_occ),
                                "_cut_",cut_label,"",
                                ".rds")    
    if(flag_redo_level>=8){
      if(cut_label=="random") {
        set.seed(1)
        fold <- dismo::kfold(myocc_sp, k=rep_all)
        myocc_sp$cut= fold
        cut_bumber <- max(fold) 
      }
      if(cut_label=="cluster") {
        set.seed(1)
        eb <- blockCV::envBlock(rasterLayer = wc_raw[[1]],
                                speciesData = myocc_sp,
                                species = NULL, 
                                k = rep_all,
                                standardization = "normal",
                                rasterBlock = FALSE,
                                numLimit = 50)        
        myocc_sp$cut= eb$foldID
        cut_bumber <- max(eb$foldID) 
        if(F){
          ggplot(myocc_sp@data,aes(bio10,bio17,col=factor(cut) ))+geom_point()
        }
      }
      if(cut_label=="block"){
        cut <- my_blockcut (   occ=as.data.frame(myocc_sp@coords)   ) 
        myocc_sp$cut= cut$occ.grp
        cut_bumber <- max(cut$bg.grp)
      }
      if(cut_label=="block2"){
        sb <- blockCV::spatialBlock(
          speciesData = myocc_sp,
          species = NULL, 
          k = rep_all,
          rows = rep_all,
          cols=1,          
          numLimit = 50)        
        myocc_sp$cut= sb$foldID
        cut_bumber <- max(sb$foldID)
      }
      if(F){
        plot(myocc_sp)
        plot(subset(myocc_sp,cut==1),col=1,add=T )
        plot(subset(myocc_sp,cut==2),col=2,add=T )
        plot(subset(myocc_sp,cut==3),col=3,add=T )
        plot(subset(myocc_sp,cut==4),col=4,add=T )
      }
      saveRDS(myocc_sp, file=path_temp_occ_cut)
    } else {
      cat("load occ; ")
      myocc_sp <- readRDS(path_temp_occ_cut)
    }    
    rep_label_comb <- combn(1:rep_all, 1)
    rep_label=1
    cl<-makeCluster(cup_num)
    registerDoParallel(cl)
    foreach(rep_label = 1:ncol(rep_label_comb),
            .combine=rbind,
            .packages=c("raster","dismo","pROC","ENMeval","maxnet"),
            .export=c("catf","calPAUC","get.params","calc.aicc","myENM_auc")
    ) %dopar% {
      rep_sel <- rep_label_comb[,rep_label]
      occ_train <- myocc_sp[which( !myocc_sp$cut %in% rep_sel),]
      if(flag_testP=="testOnly"){
        occ_test <-  myocc_sp[which(   myocc_sp$cut  %in% rep_sel),] 
      } else {
        occ_test <-  myocc_sp
      }      
      occ_train <- occ_train[names(wc_raw)]
      occ_test <- occ_test[names(wc_raw)]
      buff_label=5000
      for(buff_label in buffersize_all){
        print(buff_label)
        if(T){
          path_temp_buffer <- paste0(path_input_occ_buff,"/",
                                     basename(this_species$csv_occ),
                                     "_cut_",cut_label,"",
                                     "_rep_",rep_label,"",
                                     "_buff_",buff_label,"",
                                     ".rds")
          path_temp_buffer_env <- paste0(path_input_occ_buff,"/",
                                         basename(this_species$csv_occ),
                                         "_cut_",cut_label,"",
                                         "_rep_",rep_label,"",
                                         "_buff_",buff_label,"_env",
                                         "bgnum",GLOBAL_bgnum,
                                         ".rds")
          path_temp_buffer_env_layers <- paste0(path_input_occ_buff,"/",
                                                basename(this_species$csv_occ),
                                                "_cut_",cut_label,"",
                                                "_rep_",rep_label,"",
                                                "_buff_",buff_label,"_envLayers",
                                                "bgnum",GLOBAL_bgnum,
                                                
                                                ".rds")
          path_temp_map <- paste0(path_input_map,"/",
                                  basename(this_species$csv_occ),
                                  "_cut_",cut_label,"",
                                  "_rep_",rep_label,"",
                                  "_buff_",buff_label,
                                  "bgnum",GLOBAL_bgnum,
                                  ".tiff")
          path_temp_model <- paste0(path_output_model,"/",
                                    basename(this_species$csv_occ),
                                    "_cut_",cut_label,"",
                                    "_rep_",rep_label,"",
                                    "_buff_",buff_label,
                                    "_feature_",MYFEATURES,
                                    "_clamp_",flag_clamp,
                                    "_testP_",flag_testP,
                                    "_testA_",flag_testA,
                                    "bgnum",GLOBAL_bgnum,
                                    ".rds")
          
          path_temp_eval <- paste0(path_output_eval,"/",
                                   basename(this_species$csv_occ),
                                   "_cut_",cut_label,"",
                                   "_rep_",rep_label,"",
                                   "_buff_",buff_label,
                                   "_feature_",MYFEATURES,
                                   "_clamp_",flag_clamp,
                                   "_testP_",flag_testP,
                                   "_testA_",flag_testA,
                                   "bgnum",GLOBAL_bgnum,
                                   ".rds")
        }
        if(buff_label==-9999){
          one_buff <- shapefile(this_species$path_M)
        } else{
          if(buffer_method=="allOcc"){
            one_buff <- buffer(myocc_sp,buff_label*1000 )
          } else{
            one_buff <- buffer(occ_train,buff_label*1000 )
          }
          catf(paste0("build occ","\n"))
        }
        saveRDS(one_buff,file=path_temp_buffer)
        wc_buffer <- crop(wc_raw,extent(one_buff))
        wc_buffer <- mask(wc_buffer, one_buff)
        saveRDS(wc_buffer,file = path_temp_buffer_env_layers)       
        wc_buffer_df <- as.data.frame(wc_buffer,na.rm=TRUE,xy=TRUE)
        max_cell <- nrow(wc_buffer_df)        
        if(GLOBAL_bgnum=="all"){
          one_bg <- wc_buffer_df
        } else{
          if(max_cell < GLOBAL_bgnum){
            this_bgnum <- max_cell
          } else {
            this_bgnum <- GLOBAL_bgnum
          }
          set.seed(1) 
          one_bg <- wc_buffer_df[sample(1:nrow(wc_buffer_df),
                                        size = this_bgnum),]
        }
        saveRDS(one_bg,file = path_temp_buffer_env)  
        if(flag_testA=="bioM"){
          bg_trainMAXBG <- readRDS(gsub("_buff_[0-9]+_", 
                                        paste0("_buff_",-9999,"_"),
                                        path_temp_buffer_env))
        } else {
          bg_trainMAXBG <- readRDS(gsub("_buff_[-]?[0-9]+_",
                                        paste0("_buff_",MAXBUFFER,"_"),
                                        path_temp_buffer_env))          
        }
        temp_r <- !is.na(wc_buffer[[1]])
        maskbg_cell_num <- sum(values(temp_r))
        one_buffer_size <- raster::area(one_buff)/1000000         
        tiff(path_temp_map)
        plot(wc_buffer[[1]])
        plot(one_buff,add=T)
        points(one_bg$x[1:1000],one_bg$y[1:1000])
        plot(occ_train,add=T,col="red")
        plot(occ_test,add=T,col="blue")
        dev.off()
        bg_trainMAXBG <-  bg_trainMAXBG[! names(bg_trainMAXBG) %in% c("x","y")]
        bg_train <-  one_bg[! names(one_bg) %in% c("x","y")]
        ap <- c(rep(1,nrow(occ_train)),rep(0,nrow(bg_train)) )
        pder <- rbind(occ_train@data,bg_train)
        for(which_algorothm in all_algorithm){
          if(which_algorothm=="glm"){
            p.wt = rep(1.e-6, length(ap)) 
            p.wt[ap == 0] = 10000/sum(ap == 0) 
            
            df_glm <- cbind(ap,pder)
            one_model <- glm(ap~bio10+bio11+bio16+bio17 +
                               I(bio10^2)+I(bio11^2)+I(bio16^2)+I(bio17^2),
                             family=poisson() ,data=df_glm,
                             weights = p.wt)            
          }
          if(which_algorothm=="maxnet"){
            tryCatch({
              one_model <- maxnet::maxnet(p=ap,data=pder,
                                          maxnet.formula( p=ap,data=pder,
                                                          classes=MYFEATURES)  )
            }, error=function(e){cat("ERROR :",
                                     conditionMessage(e),
                                     "\n")})
            
          }
          save(one_model,file=path_temp_model)         
          one_eval       <- myENM_auc(one_model,
                                      occ_train@data,
                                      occ_test@data,
                                      bg_train=bg_train,
                                      bg_test=bg_train,
                                      flag_clamp=flag_clamp,
                                      which_algorothm = which_algorothm)
          one_eval_MAXBG <- myENM_auc(one_model,
                                      occ_train@data, 
                                      occ_test@data,
                                      bg_train=bg_train,
                                      bg_test=bg_trainMAXBG,
                                      flag_clamp=flag_clamp,
                                      which_algorothm = which_algorothm)
          colnames(one_eval_MAXBG) <- paste0("MAXBG",colnames(one_eval_MAXBG))
          one_line <- cbind(this_species,
                            "algorothm"=which_algorothm,
                            "cut"=cut_label,
                            "path_temp_occ_cut"=path_temp_occ_cut,
                            "path_temp_buffer"=path_temp_buffer,
                            "path_temp_buffer_env"=path_temp_buffer_env,
                            path_temp_map=path_temp_map,
                            path_temp_model=path_temp_model,
                            "fold"=rep_label,
                            "buff"=buff_label,
                            "bufferarea"=one_buffer_size,
                            "bufferarea_ncell"=maskbg_cell_num,
                            "MYFEATURES"=MYFEATURES,
                            "buffer_method"=buffer_method,
                            one_eval,
                            one_eval_MAXBG)
          saveRDS(one_line,path_temp_eval)
        }
      }
    }
    stopCluster(cl)
    gc()
  }
} 

library(raster)
library(dismo)
library(blockCV)
library(maxnet)
library(foreach)
library(doParallel)
path_input_index <- "a_runFolder/index"
occ_summary <- readRDS(paste0(path_input_index,"/occ_summary_v1.rds"))
for(mm in 1:nrow(occ_summary)  ){
    myWorkflow(occ_summary[mm,], cup_num=15,
               flag_clamp=FALSE,
               cut_all=c("cluster","random"), 
               path_input = "b_runFolder" ,
               buffer_method = "allOcc",
               flag_testA = "largest",
               buffersize_all =c(seq(5000,1000,-1000),
                                 seq(900,100,-100),
                                 seq(90,10,-10),
                                 seq(9,5,-2),
                                 -9999)
    )
  }
for(mm in 1:nrow(occ_summary)  ){
    myWorkflow(occ_summary[mm,], cup_num=15,
               flag_clamp=FALSE,
               cut_all=c("cluster","random"),  
               path_input = "b_runFolder_bioM" ,
               buffer_method = "allOcc",
               flag_testA = "bioM",
               buffersize_all =c(-9999,
                                 seq(5000,1000,-1000),
                                 seq(900,100,-100),
                                 seq(90,10,-10),
                                 seq(9,5,-2)
               ))
  }


