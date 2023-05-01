
library(ggplot2)
library(data.table)
path_experiment = c("b_runFolder_bioM")
path_experiment = c("b_runFolder")
path_output_fig = paste0(path_experiment,"/output_fig")
dir.create(path_output_fig[1])
path_output_eval = paste0(path_experiment,"/output_eval")
files_eval <- list.files(path_output_eval[1],full.names = T)
files_eval1 <- list.files(path_output_eval[2],full.names = T)
files_eval = c(files_eval,files_eval1)
length(files_eval)
all_eval <- vector()
for(i in 1:length(files_eval) ){
  print(i)
  one_record <- readRDS(files_eval[i])
  if(!"n_train_occ" %in%   names(one_record)){
    one_record$n_train_occ=-9999
    one_record$n_test_occ=-9999
    one_record$MAXBGn_train_occ=-9999
    one_record$MAXBGn_test_occ=-9999
    one_record <- one_record[names(all_eval)]
  } 
  if( grepl("clamp_FALSE",files_eval[i]) ){
    one_record$flag_clamp <- "NoClamp"
  } else {
    one_record$flag_clamp <- "YesClamp"
  }
  if( grepl("bgnum10000",files_eval[i]) ){
    one_record$trainingBG <- "10000"
  } else {
    one_record$trainingBG <- "all"
  }
  all_eval <- rbind(all_eval,one_record)
}
print(nrow(all_eval))
  df <- as.data.frame(all_eval)
  setDT(df)
  df$areatype = "buffer"
  df$areatype[df$buff==-9999] <- "bioM"
  sss=1
  for(sss in 1: length( unique(df$species) ) ){
    this_spp <- unique(df$species)[sss]
    dd_spp <- subset(df,species==this_spp)
    are_diff <- abs( dd_spp$bufferarea_ncell - 
                       unique(  dd_spp$bufferarea_ncell[dd_spp$buff=="-9999"])   )
    are_diff[are_diff==0] <- max(are_diff)
    select_size <- which.min(are_diff)
    size_most_close <- dd_spp$buff[select_size]
    df$buff[df$species==this_spp &  df$buff==-9999] <- size_most_close
  }
  df$flag_color <- paste(df$cut,df$areatype,sep="_")
saveRDS(df,paste0(path_output_fig,"/merged_df.rds"))

library(ggplot2)
library(data.table)
library(segmented)
path_experiment = c("b_runFolder_bioM")
path_experiment = c("b_runFolder")
path_output_fig = paste0(path_experiment,"/output_fig")
dir.create(path_output_fig[1])
df = readRDS(paste0(path_output_fig,"/merged_df.rds"))
flag_rescale = T
flag_rescale = F
if(flag_rescale){ 
  scale01 = function(x){  (x-mean(x))/sd(x)      }
  df_copy = copy(df)  
  df[, `:=` (
    sensitivity5 = scale01(sensitivity5),
    MAXBGspecificity5 = scale01(MAXBGspecificity5),
    MAXBGTSS5 = scale01(MAXBGTSS5),
    MAXBGboyce_test = scale01(MAXBGboyce_test)  ), 
    by = species]
}
path_output_regression = paste0(path_experiment,"/output_regression")
dir.create(path_output_regression[1])
p1= Ttest_all("sensitivity5",paired_Ttest=T)
p2= Ttest_all("MAXBGspecificity5",paired_Ttest=T)
p3= Ttest_all("MAXBGTSS5",paired_Ttest=T)
p4= Ttest_all("MAXBGboyce_test",paired_Ttest=T)
out_template = data.frame(
  row.names=c("sensitivity5","sensitivity5_greater","sensitivity5_less",
              "MAXBGspecificity5","MAXBGspecificity5_greater","MAXBGspecificity5_less",
              "MAXBGTSS5", "MAXBGTSS5_greater", "MAXBGTSS5_less",
              "MAXBGboyce_test","MAXBGboyce_test_greater","MAXBGboyce_test_less"),
           random_table2=rep(NA,12),env_table2=rep(NA,12),
           random_table3=rep(NA,12),env_table3=rep(NA,12)
  )
out_Table = out_template
out_Table[1:3,1:4] = matrix( unlist( p1[1,2:13]), nrow=3, ncol=4)
out_Table[4:6,1:4] = matrix( unlist( p2[1,2:13]), nrow=3, ncol=4)
out_Table[7:9,1:4] = matrix( unlist( p3[1,2:13]), nrow=3, ncol=4)
out_Table[10:12,1:4] = matrix( unlist( p4[1,2:13]), nrow=3, ncol=4)
write.csv(out_Table,paste0(path_output_regression,"/t_test.csv")   )
temp_scatter_df = readRDS(paste0(path_output_regression[1],
                                         "/turning_point","_",
                                         "allSPP","_",
                                         "2cut","_",
                                         "4var","_",
                                         "flag_rescale_",flag_rescale,
                                         ".rds") )
Ttest_all = function(this_var="sensitivity5",paired_Ttest=FALSE){
  df_Ttest=data.frame(species=rep(NA,length(unique(df$species) )*2 ),
                      cut=NA,
                      bioM_mean=NA,
                      buffer_mean=NA,
                      p_value=NA,
                      p_value_greater=NA,
                      buff=NA)
  for(this_sp in unique(df$species)){
    for(this_cut in unique(df$cut)){
      df_compare1 = df[species==this_sp & cut==this_cut]      
      size_to_compare =temp_scatter_df[species==this_sp &  cut==this_cut & index ==this_var]$turningP
      g2 = df_compare1[buff>=size_to_compare & areatype!= "bioM",  ..this_var ]
      g1 = df_compare1[areatype=="bioM",..this_var]
      t_results = t.test(unlist(g1),unlist(g2),paired=paired_Ttest)
      t_results_greater = t.test(unlist(g1),unlist(g2),alternative = "greater",paired=paired_Ttest)
      t_results_less = t.test(unlist(g1),unlist(g2),alternative = "less",paired=paired_Ttest)
      df_Ttest$species[i]    =this_sp
      df_Ttest$cut[i]        =this_cut
      df_Ttest$bioM_mean[i]  =t_results$estimate[1]
      df_Ttest$buffer_mean[i]=t_results$estimate[2]
      df_Ttest$p_value[i]    =t_results$p.value
      df_Ttest$p_value_greater[i]    =t_results_greater$p.value
      df_Ttest$p_value_less[i]    =t_results_less$p.value
      if(is.na(df_Ttest$p_value[i] )){
        df_Ttest$p_value[i]=999
      }
      if(is.na(df_Ttest$p_value_greater[i] )){
        df_Ttest$p_value_greater[i]=999
      }
      if(is.na(df_Ttest$p_value_less[i] )){
        df_Ttest$p_value_less[i]=999
      }
      df_Ttest$buff[i]    =size_to_compare
      i=i+1
    }
  }
  return(
    data.frame(this_var=this_var,
               n_random_t = sum(df_Ttest[df_Ttest$cut=="random" ,]$p_value<0.05 ),
               n_random_tgreater = sum(df_Ttest[df_Ttest$cut=="random" ,]$p_value_greater<0.05 ),
               n_random_tless = sum(df_Ttest[df_Ttest$cut=="random" ,]$p_value_less<0.05 ),
               n_env_t = sum(df_Ttest[df_Ttest$cut=="cluster" ,]$p_value<0.05 ),
               n_env_tgreater = sum(df_Ttest[df_Ttest$cut=="cluster",]$p_value_greater<0.05 ),
               n_env_tless = sum(df_Ttest[df_Ttest$cut=="cluster" ,]$p_value_less<0.05 )
    )
  )
}
p1= Ttest_all("sensitivity5",paired_Ttest=F)
p2= Ttest_all("MAXBGspecificity5",paired_Ttest=F)
p3= Ttest_all("MAXBGTSS5",paired_Ttest=F)
p4= Ttest_all("MAXBGboyce_test",paired_Ttest=F)

out_template = data.frame(
  row.names=c("sensitivity5","sensitivity5_greater","sensitivity5_less",
              "MAXBGspecificity5","MAXBGspecificity5_greater","MAXBGspecificity5_less",
              "MAXBGTSS5", "MAXBGTSS5_greater", "MAXBGTSS5_less",
              "MAXBGboyce_test","MAXBGboyce_test_greater","MAXBGboyce_test_less"),
  random_table2=rep(NA,12),env_table2=rep(NA,12),
  random_table3=rep(NA,12),env_table3=rep(NA,12)
)
out_Table = out_template
out_Table[1:3,1:4] = matrix( unlist( p1[1,2:7]), nrow=3, ncol=4)
out_Table[4:6,1:4] = matrix( unlist( p2[1,2:7]), nrow=3, ncol=4)
out_Table[7:9,1:4] = matrix( unlist( p3[1,2:7]), nrow=3, ncol=4)
out_Table[10:12,1:4] = matrix( unlist( p4[1,2:7]), nrow=3, ncol=4)
write.csv(out_Table,   paste0(path_output_regression,"/t_test_bioM_vs_largerBuffers_debug.csv")   )

model_list = list(); i=1
tag_list = vector()
for(sel_var in c("sensitivity5","MAXBGspecificity5","MAXBGTSS5","MAXBGboyce_test")  ){
    for(sel_cut in unique(df$cut)){
    print(i)
    sel_columns = c("species",sel_var,"buff","fold","areatype")
    part1 = df[cut==sel_cut& areatype=="bioM",..sel_columns]
    part2 = df[cut==sel_cut& areatype!="bioM",..sel_columns]    
    part12= merge(part1,part2,by=c("species","fold","buff" ),all.x=TRUE)   
    unsel_columns = paste0(c(sel_var,"areatype"),".x")
    part12[,(unsel_columns):=NULL]
    part2_new = part12
    names(part2_new)  = gsub("[.]y" , "",names(part2_new))
    setcolorder(part2_new,names(part1))
    part12_new  = rbind(part1,part2_new)
    library(lme4)
    library(lmerTest)    
    formula=   as.formula( paste( sel_var, "~", 
                                  paste("areatype", collapse = "+"),
                                  "+(1|species)"))
    model_list[[i]] = lmer(formula, data = part12_new)     
    tag_list[i] = paste(sel_var,sel_cut)
    i=i+1
  }
}

myTable <- texreg::htmlreg(l= model_list, custom.model.names=tag_list,	
                   digits=2,bold=0.05,
                   caption="Table X. Model summary",caption.above=TRUE,
                   include.aic = TRUE,single.row=FALSE,
                   file=paste0(path_output_regression,"/mixed_model_bioM_vs_similarSize",
                               "_flag_rescale_",flag_rescale,"_scaleto0mean1sd",
                               ".html"))

model_list = list(); i=1
tag_list = vector()
for(sel_var in c("sensitivity5","MAXBGspecificity5","MAXBGTSS5","MAXBGboyce_test")  ){
  for(sel_cut in unique(df$cut)){
    print(i)
    sel_columns = c("species",sel_var,"buff","fold","areatype")
    part1 = df[cut==sel_cut& areatype=="bioM",..sel_columns]
    part2 = df[cut==sel_cut& areatype!="bioM",..sel_columns]
    part12= merge(part1,part2,by=c("species","fold" ),all.x=TRUE)
    part12= part12[buff.y>=buff.x,]
    part12 = part12[,-c("buff.y")]
    names(part12)[4] = "buff"
    unsel_columns = paste0(c(sel_var,"areatype"),".x")
    part12[,(unsel_columns):=NULL]
    part2_new = part12
    names(part2_new)  = gsub("[.]y" , "",names(part2_new))
    setcolorder(part2_new,names(part1))
    part12_new  = rbind(part1,part2_new)
    library(lme4)
    library(lmerTest)    
    formula=   as.formula( paste( sel_var, "~", 
                                  paste("areatype", collapse = "+"),
                                  "+(1|species)"))
    model_list[[i]] = lmer(formula, data = part12_new)     
    tag_list[i] = paste(sel_var,sel_cut)
    i=i+1
  }
}
myTable <- texreg::htmlreg(l= model_list, custom.model.names=tag_list,	
                           digits=2,bold=0.05,
                           caption="Table X. Model summary",caption.above=TRUE,
                           include.aic = TRUE,single.row=FALSE,
                           file=paste0(path_output_regression,"/mixed_model_bioM_vs_biggerSize",
                                       "_flag_rescale_",flag_rescale,"_scaleto0mean1sd",
                                       ".html"))
temp_scatter_df = readRDS(paste0(path_output_regression[1],
                                 "/turning_point","_",
                                 "allSPP","_",
                                 "2cut","_",
                                 "4var","_",
                                 "flag_rescale_",flag_rescale,
                                 ".rds") )
model_list = list(); i=1
tag_list = vector()
for(sel_var in c("sensitivity5","MAXBGspecificity5","MAXBGTSS5","MAXBGboyce_test")  ){
  for(sel_cut in unique(df$cut)){
    print(i)
    sel_columns = c("species",sel_var,"buff","fold","areatype")
    part1 = df[cut==sel_cut& areatype=="bioM",..sel_columns]
    part2 = df[cut==sel_cut& areatype!="bioM",..sel_columns]
    temp_scatter_df_subset = temp_scatter_df[cut == sel_cut & index == sel_var]
    part2_temp = merge(part2,temp_scatter_df_subset[,c("species","turningP")],by="species", all.x=T)
    part2 =part2_temp[buff>=turningP]
    part2 = part2[,-"turningP"]
    part12= merge(part1,part2,by=c("species","fold" ),all.x=TRUE)
    part12= part12[buff.y>=buff.x,]
    part12 = part12[,-c("buff.y")]
    names(part12)[4] = "buff"
    unsel_columns = paste0(c(sel_var,"areatype"),".x")
    part12[,(unsel_columns):=NULL]
    part2_new = part12
    names(part2_new)  = gsub("[.]y" , "",names(part2_new))
    setcolorder(part2_new,names(part1))
    part12_new  = rbind(part1,part2_new)
    library(lme4)
    library(lmerTest)    
    formula=   as.formula( paste( sel_var, "~", 
                                  paste("areatype", collapse = "+"),
                                  "+(1|species)"))
    model_list[[i]] = lmer(formula, data = part12_new)     
    tag_list[i] = paste(sel_var,sel_cut)
    i=i+1
  }
}
myTable <- texreg::htmlreg(l= model_list, custom.model.names=tag_list,	
                           digits=2,bold=0.05,
                           caption="Table X. Model summary",caption.above=TRUE,
                           include.aic = TRUE,single.row=FALSE,
                           file=paste0(path_output_regression,"/mixed_model_bioM_vs_2segment",
                                       "_flag_rescale_",flag_rescale,"_scaleto0mean1sd",
                                       ".html"))
model_list = list(); i=1
tag_list = vector()
sel_var="sensitivity5";
for(sel_var in c("sensitivity5","MAXBGspecificity5","MAXBGTSS5","MAXBGboyce_test")  ){
  for(sel_cut in unique(df$cut)){
    print(i)
    sel_columns = c("species",sel_var,"buff","fold","areatype")
    part1 = df[cut==sel_cut& areatype=="bioM",..sel_columns]
    part2 = df[cut==sel_cut& areatype!="bioM",..sel_columns]
    df_compare_temp <-   part2[,c(  lapply(.SD,mean)) ,
                                     by=c("areatype","buff","species"),
                                     .SDcols=c(  sel_var   )]
    names(df_compare_temp)[4] = "sel_index"
    setorder(df_compare_temp, species, -sel_index)
    df_compare_temp <- unique(df_compare_temp, by = "species")
    part2_updated = merge(part2,df_compare_temp,by=c("areatype","buff","species"), all.y=T)
    part2_updated = part2_updated [,-"sel_index"]
    setcolorder(part2_updated,names(part1))
    part12_new  = rbind(part1,part2_updated)
    library(lme4)
    library(lmerTest)
    formula=   as.formula( paste( sel_var, "~", 
                                  paste("areatype", collapse = "+"),
                                  "+(1|species)"))
    model_list[[i]] = lmer(formula, data = part12_new) 
    tag_list[i] = paste(sel_var,sel_cut)
    i=i+1
  }
}
myTable <- texreg::htmlreg(l= model_list, custom.model.names=tag_list,	
                           digits=2,bold=0.05,
                           caption="Table X. Model summary",caption.above=TRUE,
                           include.aic = TRUE,single.row=FALSE,
                           file=paste0(path_output_regression,"/mixed_model_bioM_vs_best",
                                       "_flag_rescale_",flag_rescale,"_scaleto0mean1sd",
                                       ".html"))

this_var <- "sensitivity5" # "npap" #
testdf <-   df[,c(  lapply(.SD,mean),
                    lapply(.SD,sd)
) ,
by=c("areatype","cut","buff","MYFEATURES","buffer_method","species","flag_clamp" ),
.SDcols=c(  this_var   )]
name_index = which(names(testdf)==this_var)
names(testdf)[name_index[1]] <- "dependent_var"
names(testdf)[name_index[2]]  <- "dependent_var_sd"
this_sp = unique(df$species)[2]
this_sp = unique(df$species)[13]
for(this_sp in unique(df$species)){
  for(this_cut in unique(df$cut)){
    this_df =testdf[cut==this_cut & 
                      species==this_sp & 
                      buff <=2000 &
                      areatype=="buffer"]
    bioM_df =testdf[cut==this_cut & 
                      species==this_sp & 
                      buff <=2000 &
                      areatype=="bioM"]    
    plot(this_df$buff,this_df$dependent_var)    
    flag_log10=TRUE
    if(flag_log10){
    this_df$buff = log10(this_df$buff)
    guess_x=2
    plot(this_df$buff,this_df$dependent_var)
    } else{
      guess_x=100
    }
    this_m <- lm(dependent_var ~ buff, data = this_df)
    o <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                   control = seg.control(display = FALSE)    )
    plot(o)
    print(o$psi[2]) 
    dat2 = data.frame(buff = this_df$buff, dependent_var = broken.line(o)$fit)
    plot(dat2$buff,dat2$dependent_var)    
    if(flag_log10){
      this_df$buff = 10^(this_df$buff)
      dat2$buff = 10^(dat2$buff)
    }    
    fig_segment = ggplot(this_df, aes(x = buff, y = dependent_var)) +
      geom_point(col="gray30") +
      geom_line(data = dat2, color = 'black')+
      geom_point(data=bioM_df,col="red")
    fig_segment
    ggsave(plot = fig_segment,
           paste0(path_output_regression[1],
                  "/segment","_",
                  this_sp,"_",
                  this_cut,"_",
                  this_var,"_",
                  ".tiff"),
           width = 10,height = 10,units = "cm",dpi = 300)    
  }
}

prep_data_2023 = function( this_cut = "random",
                      this_var = c("sensitivity5","MAXBGspecificity5","MAXBGTSS5","MAXBGboyce_test")
){
  testdf <-   df[cut==this_cut & areatype=="buffer",
                 c(  lapply(.SD,mean),
                     lapply(.SD,sd)) ,
                 by=c("areatype","cut","buff","MYFEATURES","buffer_method","flag_clamp" ),
                 .SDcols=c(  this_var   )]
  names(testdf)[11:14] = paste0(names(testdf)[11:14],"_sd")
  this_df =testdf
  testdf <-   df[cut==this_cut & areatype=="bioM",
                 c(  lapply(.SD,mean),
                     lapply(.SD,sd)) ,
                 by=c("areatype","cut","MYFEATURES","buffer_method","flag_clamp" ),
                 .SDcols=c(  this_var   )]
  names(testdf)[10:13] = paste0(names(testdf)[10:13],"_sd")
  bioM_df =testdf  
  flag_log10=TRUE
  if(flag_log10){
    this_df$buff = log10(this_df$buff)
    guess_x=2
    plot(this_df$buff,this_df$dependent_var)
  } else{
    guess_x=100
  }   
  this_m <- lm(sensitivity5 ~ buff, data = this_df)
  o1 <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                 control = seg.control(display = FALSE)    )
  dat2_sensitivity5 = data.frame(buff = this_df$buff, dependent_var = broken.line(o1)$fit)
  seg_pred_value = predict(o1,newdata=data.frame(buff=o1$psi[2]))
  dat2_sensitivity5 = rbind(dat2_sensitivity5,data.frame(buff =o1$psi[2],
                                                         dependent_var=seg_pred_value))
  this_m <- lm(MAXBGspecificity5 ~ buff, data = this_df)
  o2 <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                 control = seg.control(display = FALSE)    )
  dat2_MAXBGspecificity5 = data.frame(buff = this_df$buff, 
                                      dependent_var = broken.line(o2)$fit)
  seg_pred_value = predict(o2,newdata=data.frame(buff=o2$psi[2]))
  dat2_MAXBGspecificity5 = rbind(dat2_MAXBGspecificity5,data.frame(buff =o2$psi[2],
                                                                   dependent_var=seg_pred_value))
  this_m <- lm(MAXBGTSS5 ~ buff, data = this_df)
  o3 <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                 control = seg.control(display = FALSE)    )
  dat2_MAXBGTSS5 = data.frame(buff = this_df$buff, dependent_var = broken.line(o3)$fit)
  seg_pred_value = predict(o3,newdata=data.frame(buff=o3$psi[2]))
  dat2_MAXBGTSS5 = rbind(dat2_MAXBGTSS5,data.frame(buff =o3$psi[2],
                                                   dependent_var=seg_pred_value))
  this_m <- lm(MAXBGboyce_test ~ buff, data = this_df)
  o4 <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                 control = seg.control(display = FALSE)    )
  dat2_MAXBGboyce = data.frame(buff = this_df$buff, dependent_var = broken.line(o4)$fit)
  seg_pred_value = predict(o4,newdata=data.frame(buff=o4$psi[2]))
  dat2_MAXBGboyce = rbind(dat2_MAXBGboyce,data.frame(buff =o4$psi[2],
                                                   dependent_var=seg_pred_value))
  if(flag_log10){
    this_df$buff = 10^(this_df$buff)
    dat2_sensitivity5$buff = 10^(dat2_sensitivity5$buff)
    dat2_MAXBGspecificity5$buff = 10^(dat2_MAXBGspecificity5$buff)
    dat2_MAXBGTSS5$buff = 10^(dat2_MAXBGTSS5$buff)
    dat2_MAXBGboyce$buff = 10^(dat2_MAXBGboyce$buff)
  }  
  return(list(this_df,bioM_df,
              dat2_sensitivity5,dat2_MAXBGspecificity5,dat2_MAXBGTSS5,dat2_MAXBGboyce,
              list(o1,o2,o3,o4)  )   )
}
df_random  = prep_data_2023 (this_cut = "random")
df_cluster = prep_data_2023 (this_cut = "cluster")
buffer_df = rbind(df_random[[1]],df_cluster[[1]])
bioM_df = rbind(df_random[[2]],df_cluster[[2]])
df_random[[3]]$cut="random"
df_cluster[[3]]$cut="cluster"
segment_df_sen = rbind(df_random[[3]],df_cluster[[3]])
df_random[[4]]$cut="random"
df_cluster[[4]]$cut="cluster"
segment_df_spec = rbind(df_random[[4]],df_cluster[[4]])
df_random[[5]]$cut="random"
df_cluster[[5]]$cut="cluster"
segment_df_tss = rbind(df_random[[5]],df_cluster[[5]])
df_random[[6]]$cut="random"
df_cluster[[6]]$cut="cluster"
segment_df_boyce = rbind(df_random[[6]],df_cluster[[6]])
df_random[[3]]$buff[nrow(df_random[[3]])]
df_cluster[[3]]$buff[nrow(df_cluster[[3]])]
df_random[[4]]$buff[nrow(df_random[[4]])]
df_cluster[[4]]$buff[nrow(df_cluster[[4]])]
df_random[[5]]$buff[nrow(df_random[[5]])]
df_cluster[[5]]$buff[nrow(df_cluster[[5]])]
df_random[[6]]$buff[nrow(df_random[[6]])]
df_cluster[[6]]$buff[nrow(df_cluster[[6]])]
seg_slope_output = rbind( slope(df_random[[7]][[1]])$buff[,1],
slope(df_random[[7]][[2]])$buff[,1],
slope(df_random[[7]][[3]])$buff[,1],
slope(df_random[[7]][[4]])$buff[,1])
seg_slope_output = data.frame(seg_slope_output)
seg_slope_output$turning_p = c(
  df_random[[3]]$buff[nrow(df_random[[3]])],
  df_random[[4]]$buff[nrow(df_random[[4]])],
  df_random[[5]]$buff[nrow(df_random[[5]])],
  df_random[[6]]$buff[nrow(df_random[[6]])]   )
seg_slope_output$r2 =  c( summary(df_random[[7]][[1]])$r.squared,
summary(df_random[[7]][[2]])$r.squared,
summary(df_random[[7]][[3]])$r.squared,
summary(df_random[[7]][[4]])$r.squared)
seg_slope_output2 = rbind( slope(df_cluster[[7]][[1]])$buff[,1],
                          slope(df_cluster[[7]][[2]])$buff[,1],
                          slope(df_cluster[[7]][[3]])$buff[,1],
                          slope(df_cluster[[7]][[4]])$buff[,1])
seg_slope_output2 = data.frame(seg_slope_output2)
seg_slope_output2$turning_p = c(
  df_cluster[[3]]$buff[nrow(df_cluster[[3]])],
  df_cluster[[4]]$buff[nrow(df_cluster[[4]])],
  df_cluster[[5]]$buff[nrow(df_cluster[[5]])],
  df_cluster[[6]]$buff[nrow(df_cluster[[6]])]   )
seg_slope_output2$r2 =  c( summary(df_cluster[[7]][[1]])$r.squared,
                          summary(df_cluster[[7]][[2]])$r.squared,
                          summary(df_cluster[[7]][[3]])$r.squared,
                          summary(df_cluster[[7]][[4]])$r.squared)
seg_slope_output = cbind(seg_slope_output,seg_slope_output2)
write.csv(seg_slope_output,row.names = F,
          paste0(path_output_regression[1],
                 "/slopes","_",
                 "allSPP","_",
                 "2cut","_",
                 "4var","_",
                 "flag_rescale_",flag_rescale,
                 ".csv")
          )
fig_sensitivity5 = ggplot(buffer_df, aes(x = buff, y = sensitivity5,col=cut)) +
  geom_point(size=1) +
  geom_line(data = segment_df_sen, aes(x = buff, y = dependent_var,col=cut)) +
  geom_hline(data=bioM_df,aes(yintercept =sensitivity5,col=cut ),linetype =2)+
  scale_x_log10(breaks=c(10,100,1000,5000))+ylim(-0.05,1)+
  scale_color_discrete(
                       labels=c("envBlock","random"))+
  xlab("Buffer size (km)")+ ylab("Sensitivity")+
  annotate("text",label="bioM",x=10,y=bioM_df$sensitivity5)+
  theme(legend.position = c(0.8, 0.23),
        legend.title = element_blank()
        )
fig_sensitivity5

fig_sensitivity5_eb = fig_sensitivity5+
  geom_errorbar(aes(ymin=sensitivity5-sensitivity5_sd,
                    ymax=sensitivity5+sensitivity5_sd),width=.1,alpha=0.3  )+
  ylim(-0.1,1)
fig_sensitivity5_eb


fig_spec= ggplot(buffer_df, aes(x = buff, y = MAXBGspecificity5,col=cut)) +
  geom_point(size=1) +
  geom_line(data = segment_df_spec, aes(x = buff, y = dependent_var,col=cut)) +
  geom_hline(data=bioM_df,aes(yintercept =MAXBGspecificity5,col=cut ),linetype =2)+
  scale_x_log10(breaks=c(10,100,1000,5000))+ylim(0.2,1)+
  xlab("Buffer size (km)")+ ylab("Specificity")+
  annotate("text",label="bioM",x=10,y=bioM_df$MAXBGspecificity5)+
  theme(legend.position='none')
fig_spec
fig_spec_eb = fig_spec +
  geom_errorbar(aes(ymin=MAXBGspecificity5-MAXBGspecificity5_sd,
                    ymax=MAXBGspecificity5+MAXBGspecificity5_sd),width=.1 ,alpha=0.3 )+
  ylim(0.1,1)
fig_spec_eb

fig_tss = ggplot(buffer_df, aes(x = buff, y = MAXBGTSS5 ,col=cut)) +
  geom_point(size=1) +
  geom_line(data = segment_df_tss, aes(x = buff, y = dependent_var,col=cut)) +
  geom_hline(data=bioM_df,aes(yintercept =MAXBGTSS5 ,col=cut ),linetype =2)+
  xlab("Buffer size (km)")+ ylab("TSS")+
  scale_x_log10(breaks=c(10,100,1000,5000))+ylim(-0.7,1)+
  annotate("text",label="bioM",x=10,y=bioM_df$MAXBGTSS5)+
  theme(legend.position='none')
fig_tss

fig_tss_eb = fig_tss + 
  geom_errorbar(aes(ymin=MAXBGTSS5-MAXBGTSS5_sd,
                    ymax=MAXBGTSS5+MAXBGTSS5_sd),width=.1  ,alpha=0.3)+
  ylim(-1,1)
fig_tss_eb

fig_boyce = ggplot(buffer_df, aes(x = buff, y = MAXBGboyce_test ,col=cut)) +
  geom_point(size=1) +
  geom_line(data = segment_df_boyce, aes(x = buff, y = dependent_var,col=cut)) +
  geom_hline(data=bioM_df,aes(yintercept =MAXBGboyce_test ,col=cut ),linetype =2)+
  xlab("Buffer size (km)")+ ylab("Boyce")+
  scale_x_log10(breaks=c(10,100,1000,5000))+ylim(-0.7,1)+
  annotate("text",label="bioM",x=10,y=bioM_df$MAXBGboyce_test)+
  theme(legend.position='none')
fig_boyce

fig_boyce_eb = fig_boyce + 
  geom_errorbar(aes(ymin=MAXBGboyce_test-MAXBGboyce_test_sd,
                    ymax=MAXBGboyce_test+MAXBGboyce_test_sd),width=.1 ,alpha=0.3 )+
  ylim(-1,1.1)
fig_boyce_eb

library(gridExtra)
combined_fig = grid.arrange(fig_sensitivity5,
                            fig_spec,
                            fig_tss,
                            fig_boyce,
                            ncol = 1)
ggsave(plot = combined_fig,
       paste0(path_output_regression[1],
              "/segment","_",
              "allSPP","_",
              "2cut","_",
              "4var","_",
              "flag_rescale_",flag_rescale,
              ".tiff"),
       width = 10,height = 20,units = "cm",dpi = 300)


combined_fig_eb = grid.arrange(fig_sensitivity5_eb,
                            fig_spec_eb,
                            fig_tss_eb,
                            fig_boyce_eb,
                            ncol = 1)
ggsave(plot = combined_fig_eb,
       paste0(path_output_regression[1],
              "/segment","_",
              "allSPP","_",
              "2cut","_",
              "4var","_",
              "flag_rescale_",flag_rescale,
              "_withErrorBar",
              ".tiff"),
       width = 10,height = 20,units = "cm",dpi = 300)
plot_1sp = function(whichsp = "",flag_plot=TRUE){
  prep_data = function( this_cut = "random",
                        this_var = c("sensitivity5",
                                     "MAXBGspecificity5",
                                     "MAXBGTSS5",
                                     "MAXBGboyce_test")
                        ){    
    testdf <-   df[cut==this_cut & areatype=="buffer" & species==whichsp,
                   c(  lapply(.SD,mean),
                       lapply(.SD,sd)) ,
                   by=c("areatype","cut","buff","MYFEATURES","buffer_method","flag_clamp" ),
                   .SDcols=c(  this_var   )]    
    names(testdf)[11:14] = paste0(names(testdf)[11:14],"_sd")
    this_df =testdf    
    testdf <-   df[cut==this_cut & areatype=="bioM"& species==whichsp,
                   c(  lapply(.SD,mean),
                       lapply(.SD,sd)) ,
                   by=c("areatype","cut","buff","MYFEATURES","buffer_method","flag_clamp" ),
                   .SDcols=c(  this_var   )]
    names(testdf)[11:14] = paste0(names(testdf)[11:14],"_sd")
    bioM_df =testdf    
    flag_log10=TRUE
    if(flag_log10){
      this_df$buff = log10(this_df$buff)
      guess_x=2
      plot(this_df$buff,this_df$dependent_var)
    } else{
      guess_x=100
    } 
    this_m <- lm(sensitivity5 ~ buff, data = this_df)
    o <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                   control = seg.control(display = FALSE)    )
    dat2_sensitivity5 = data.frame(buff = this_df$buff, dependent_var = broken.line(o)$fit)
    seg_pred_value = predict(o,newdata=data.frame(buff=o$psi[2]))
    dat2_sensitivity5 = rbind(dat2_sensitivity5,data.frame(buff =o$psi[2],
                                                           dependent_var=seg_pred_value))
    this_m <- lm(MAXBGspecificity5 ~ buff, data = this_df)
    o <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                   control = seg.control(display = FALSE)    )
    dat2_MAXBGspecificity5 = data.frame(buff = this_df$buff, dependent_var = broken.line(o)$fit)
    seg_pred_value = predict(o,newdata=data.frame(buff=o$psi[2]))
    dat2_MAXBGspecificity5 = rbind(dat2_MAXBGspecificity5,data.frame(buff =o$psi[2],
                                                                     dependent_var=seg_pred_value))
    
    this_m <- lm(MAXBGTSS5 ~ buff, data = this_df)
    o <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                   control = seg.control(display = FALSE)    )
    dat2_MAXBGTSS5 = data.frame(buff = this_df$buff, dependent_var = broken.line(o)$fit)
    seg_pred_value = predict(o,newdata=data.frame(buff=o$psi[2]))
    dat2_MAXBGTSS5 = rbind(dat2_MAXBGTSS5,data.frame(buff =o$psi[2],
                                                     dependent_var=seg_pred_value))
    
    this_m <- lm(MAXBGboyce_test ~ buff, data = this_df)
    o <- segmented(this_m, seg.Z = ~buff, psi = list(buff = c(guess_x)),
                   control = seg.control(display = FALSE)    )
    dat2_MAXBGboyce = data.frame(buff = this_df$buff, dependent_var = broken.line(o)$fit)
    seg_pred_value = predict(o,newdata=data.frame(buff=o$psi[2]))
    dat2_MAXBGboyce = rbind(dat2_MAXBGboyce,data.frame(buff =o$psi[2],
                                                     dependent_var=seg_pred_value))
    if(flag_log10){
      this_df$buff = 10^(this_df$buff)
      dat2_sensitivity5$buff = 10^(dat2_sensitivity5$buff)
      dat2_MAXBGspecificity5$buff = 10^(dat2_MAXBGspecificity5$buff)
      dat2_MAXBGTSS5$buff = 10^(dat2_MAXBGTSS5$buff)
      dat2_MAXBGboyce$buff = 10^(dat2_MAXBGboyce$buff)
    }    
    return(list(this_df,bioM_df,
                dat2_sensitivity5,dat2_MAXBGspecificity5,dat2_MAXBGTSS5,dat2_MAXBGboyce)   )
  }  
  df_random  = prep_data (this_cut = "random")
  df_cluster = prep_data (this_cut = "cluster")  
  buffer_df = rbind(df_random[[1]],df_cluster[[1]])
  bioM_df = rbind(df_random[[2]],df_cluster[[2]])  
  df_random[[3]]$cut="random"
  df_cluster[[3]]$cut="cluster"
  segment_df_sen = rbind(df_random[[3]],df_cluster[[3]])  
  df_random[[4]]$cut="random"
  df_cluster[[4]]$cut="cluster"
  segment_df_spec = rbind(df_random[[4]],df_cluster[[4]])  
  df_random[[5]]$cut="random"
  df_cluster[[5]]$cut="cluster"
  segment_df_tss = rbind(df_random[[5]],df_cluster[[5]])  
  df_random[[6]]$cut="random"
  df_cluster[[6]]$cut="cluster"
  segment_df_boyce = rbind(df_random[[6]],df_cluster[[6]])  
  df_random[[3]]$buff[nrow(df_random[[3]])]
  df_cluster[[3]]$buff[nrow(df_cluster[[3]])]
  df_random[[4]]$buff[nrow(df_random[[4]])]
  df_cluster[[4]]$buff[nrow(df_cluster[[4]])]
  df_random[[5]]$buff[nrow(df_random[[5]])]
  df_cluster[[5]]$buff[nrow(df_cluster[[5]])]
  df_random[[6]]$buff[nrow(df_random[[6]])]
  df_cluster[[6]]$buff[nrow(df_cluster[[6]])]  
  df_turningP_1sp = data.frame(
      species= whichsp,
      cut=c("random","cluster","random","cluster","random","cluster","random","cluster"),
             index=c("sensitivity5","sensitivity5",
                     "MAXBGspecificity5","MAXBGspecificity5",
                     "MAXBGTSS5","MAXBGTSS5",
                     "MAXBGboyce_test","MAXBGboyce_test"),
             turningP=c(  df_random[[3]]$buff[nrow(df_random[[3]])],
                          df_cluster[[3]]$buff[nrow(df_cluster[[3]])],
                          df_random[[4]]$buff[nrow(df_random[[4]])],
                          df_cluster[[4]]$buff[nrow(df_cluster[[4]])],
                          df_random[[5]]$buff[nrow(df_random[[5]])],
                          df_cluster[[5]]$buff[nrow(df_cluster[[5]])],
                          df_random[[6]]$buff[nrow(df_random[[6]])],
                          df_cluster[[6]]$buff[nrow(df_cluster[[6]])]
                          )
             )  
  if(flag_plot){
  fig_sensitivity5 = ggplot(buffer_df, aes(x = buff, y = sensitivity5,col=cut)) +
    geom_point(size=1) +
    geom_line(data = segment_df_sen, aes(x = buff, y = dependent_var,col=cut)) +
    geom_point(data=bioM_df,size=4,shape=8)+
    scale_x_log10(breaks=c(10,100,1000,5000))+ylim(-0.05,1)+
    scale_color_discrete(
      labels=c("envBlock","random"))+
    xlab("Buffer size (km)")+ ylab("Sensitivity")+
    annotate("text",label="bioM",x=bioM_df$buff,y=bioM_df$sensitivity5+0.05)+
    theme(legend.position = c(0.8, 0.23),
          legend.title = element_blank()
    )
  fig_sensitivity5
  
  fig_spec= ggplot(buffer_df, aes(x = buff, y = MAXBGspecificity5,col=cut)) +
    geom_point(size=1) +
    geom_line(data = segment_df_spec, aes(x = buff, y = dependent_var,col=cut)) +
    geom_point(data=bioM_df,size=4,shape=8)+
    scale_x_log10(breaks=c(10,100,1000,5000))+ylim(0.3,1)+
    xlab("Buffer size (km)")+ ylab("Specificity")+
    theme(legend.position='none')
  fig_spec  
  fig_tss = ggplot(buffer_df, aes(x = buff, y = MAXBGTSS5 ,col=cut)) +
    geom_point(size=1) +
    geom_line(data = segment_df_tss, aes(x = buff, y = dependent_var,col=cut)) +
    geom_point(data=bioM_df,size=4,shape=8)+
    xlab("Buffer size (km)")+ ylab("TSS")+
    scale_x_log10(breaks=c(10,100,1000,5000))+ylim(-0.6,1)+
    theme(legend.position='none')
  fig_tss
  
  library(gridExtra)
  combined_fig = grid.arrange(fig_sensitivity5,
                              fig_spec,
                              fig_tss, ncol = 1)
  
  ggsave(plot = combined_fig,
         paste0(path_output_regression[1],
                "/segment","_",
                whichsp,"_",
                "2cut","_",
                "3var","_",
                ".tiff"),
         width = 10,height = 20,units = "cm",dpi = 300)
  }
  return(df_turningP_1sp)
}

all_turnP_df = vector()
for(this_sp in unique(df$species)){#[-87]  
  one_df = plot_1sp(whichsp=this_sp,flag_plot=F)
  all_turnP_df=rbind(all_turnP_df,one_df)
}

head(all_turnP_df)
setDT(all_turnP_df)
this_var="turningP"
all_turnP_df[,
   c(  lapply(.SD,mean),
       lapply(.SD,sd),
       lapply(.SD,min),
       lapply(.SD,max)
       ) ,
   by=c("cut","index" ),
   .SDcols=c( this_var    )]

ggplot(all_turnP_df,aes(turningP,col=cut))+
  geom_density(n=100)+
  facet_grid(.~index)+
  scale_x_log10()

ggplot(all_turnP_df,aes(turningP,col=cut))+
  geom_boxplot()+
  facet_grid(.~index)+  scale_x_log10()

summary(all_turnP_df[cut=="random"& index=="sensitivity5"]$turningP)
summary(all_turnP_df[cut=="random"& index=="MAXBGspecificity5"]$turningP)
summary(all_turnP_df[cut=="random"& index=="MAXBGTSS5"]$turningP)
summary(all_turnP_df[cut=="random"& index=="MAXBGboyce_test"]$turningP)

summary(all_turnP_df[cut=="cluster"& index=="sensitivity5"]$turningP)
summary(all_turnP_df[cut=="cluster"& index=="MAXBGspecificity5"]$turningP)
summary(all_turnP_df[cut=="cluster"& index=="MAXBGTSS5"]$turningP)
summary(all_turnP_df[cut=="cluster"& index=="MAXBGboyce_test"]$turningP)

extracted_bioM = unique(df[areatype == "bioM",c("buff","species","areatype","cut")])
temp_scatter_df = merge(all_turnP_df,extracted_bioM,by=c("species","cut"),all.x=T)
table(temp_scatter_df$buff>=temp_scatter_df$turningP)/nrow(temp_scatter_df)

temp_scatter_df = temp_scatter_df[,-"areatype"]

saveRDS(temp_scatter_df,
        paste0(path_output_regression[1],
               "/turning_point","_",
               "allSPP","_",
               "2cut","_",
               "4var","_",
               "flag_rescale_",flag_rescale,
               ".rds") )
temp_scatter_df_melt = melt(temp_scatter_df,id.vars = c( "species","cut","index"))

temp_scatter_df_melt$index = factor(temp_scatter_df_melt$index,
                                    levels = c("sensitivity5",  
                                               "MAXBGspecificity5",
                                               "MAXBGTSS5",
                                               "MAXBGboyce_test" ) 
                                    )
levels(temp_scatter_df_melt$index) = c("Sensitivity","Specificity","TSS","Boyce")

levels(temp_scatter_df_melt$variable) = c("turning point","bioM")

temp_scatter_df_melt$cut = factor(temp_scatter_df_melt$cut,levels=c("random","cluster"))
levels(temp_scatter_df_melt$cut) = c("random" , "envBlock")

fig_turning_point = ggplot(temp_scatter_df_melt,aes(y=value,x=variable))+
  geom_boxplot(alpha=0.5) +facet_grid(index~cut)+scale_y_log10()+
  labs(y = "Size (km)",x="")

ggsave(plot = fig_turning_point,
       paste0(path_output_regression[1],
              "/turning_point","_",
              "allSPP","_",
              "2cut","_",
              "4var","_",
              "flag_rescale_",flag_rescale,
              ".tiff"),
       width = 10,height = 20,units = "cm",dpi = 300)
