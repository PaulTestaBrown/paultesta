#### Overview #####
# A set of helper functions to simulate data and estimate ACTE for binary Choice Experiments

# data_fn
# ate_fn
# acte_fn
# s_acte_fn


c_ate <- c("Control","Treatment")
c_acte_s <- c("Control","Selection")
c_acte_a <- c("Selection","Treatment")
c_select <- c("Avoid Information","Choose Information")



ates <- list(c_ate)
actes <- list(c_acte_s,c_acte_a)

ate_fn <- function(the_data, dv1="Y",d="treatment",c,weights=T,...){
  #frm<-paste(dv1,"~  age+age_sq+is_female+is_nonwhite+is_democrat+is_moderate+is_liberal+open+consc+extra+agree+neur")
  #df[,dv1] <- resid(lm(frm,df,na.action="na.exclude"))
  df <- the_data
  #print(dim(df))
  # Weights
  tmp <- as.data.frame(df[df[[d]]%in%c, ])
  #print(dim(tmp))
  if(weights==F){
    mu1 <- with(tmp, mean(tmp[tmp[[d]]==c[1], dv1],na.rm=T))
    mu2 <- with(tmp, mean(tmp[tmp[[d]]==c[2], dv1],na.rm=T))
  }
  if(weights==T){
    mu1 <- with(tmp, Hmisc::wtd.mean(tmp[tmp[[d]]==c[1], dv1],na.rm=T,weights=tmp[tmp[[d]]==c[1],"weights"]))
    mu2 <- with(tmp, Hmisc::wtd.mean(tmp[tmp[[d]]==c[2], dv1],na.rm=T,weights=tmp[tmp[[d]]==c[2],"weights"]))
  }
  diff <- mu2-mu1
  if(weights==F){
    sd1 <- with(tmp, sd(tmp[tmp[[d]]==c[1], dv1],na.rm=T))
    sd2 <- with(tmp, sd(tmp[tmp[[d]]==c[2], dv1],na.rm=T))
  }
  if(weights==T){
    sd1 <- sqrt(with(tmp, Hmisc::wtd.var(tmp[tmp[[d]]==c[1], dv1],na.rm=T,weights = tmp[tmp[[d]]==c[1],"weights"])))
    sd2 <- sqrt(with(tmp, Hmisc::wtd.var(tmp[tmp[[d]]==c[2], dv1],na.rm=T,weights = tmp[tmp[[d]]==c[2],"weights"])))
  }
  n1 <- with(tmp, sum(!is.na(tmp[tmp[[d]]==c[1], dv1])*tmp[tmp[[d]]==c[1],"weights"]))
  n2 <- with(tmp, sum(!is.na(tmp[tmp[[d]]==c[2], dv1])*tmp[tmp[[d]]==c[2],"weights"]))
  #print(c(mu1,mu2,n1,n2,sd1,sd2))
  se <- sqrt( sd1^2/n1 + sd2^2/n2)
  the_df <- (sd1^2/n1+sd2^2/n2)^2/((sd1^4)/(n1^2*(n1-1))+ (sd2^4)/(n2^2*(n2-1)))
  ll <- diff - qt(.975,the_df)*se
  ul <- diff + qt(.975,the_df)*se
  ll90 <- diff - qt(.95,the_df)*se
  ul90 <- diff + qt(.95,the_df)*se
  stat <- diff/se
  pval = 2 * pt(-abs(stat),the_df)
  result <- c(Difference = diff, SE = se, ll = ll, ul = ul,ll90=ll90,ul90=ul90, pval = pval)
  return(result)
}


# ate_fn(df, dv1="op_inds_6",d="treatment",c("Control","Selection"))
# ate_fn(df, dv1="op_inds_6",d="Treatment",c("Avoid Information","Choose Information"))



#### ACTE function ####

acte_fn <- function(dat,dv2="Y",d2="treatment",z,w,...){
  df <- dat
  N <- dim(df)[1]

  n_exp <- sum(df$C=="Experiment")
  n_choice <- sum(df$C == "Choice")
  n_ch_a <- sum(df$C=="Choice" & df$avoid01==1)



  n_select <- sum(df$C=="Choice" & df$avoid01==0, na.rm=T)
  n_avoid <- sum(df$C=="Choice" & df$avoid01==1,na.rm = T)
  #n_control <- sum(df$D_ch == "Control",na.rm = T)
  # df$weights <- rep(1,N)
  # df$weights[df$C=="Choice" & df$avoid01==0 & df$D_ch == "Control"] <-1/(n_select/n_choice)
  # df$weights[df$C=="Choice" & df$avoid01==1 & df$D_ch == "Control"] <-1/(n_control/n_avoid)
  tmp <- ate_fn(the_data=df, dv1=dv2,d = d2,c=z, weights=w)

  x <- as.numeric(tmp["Difference"])

  se_x <- as.numeric(tmp["SE"])
  if(z[1]=="Control"){
    y <- summary(lm(select01~1,df[df$C=="Choice",]))$coef[1,1]
    se_y <- summary(lm(select01~1,df[df$C=="Choice",]))$coef[1,2]}
  if(z[2]=="Treatment"){
    y <- summary(lm(avoid01~1,df[df$C=="Choice",]))$coef[1,1]
    se_y <- summary(lm(avoid01~1,df[df$C=="Choice",]))$coef[1,2]}
  mvec <- c(x=x, y= y)
  V <- diag(c(se_x,se_y)^2)
  est <- car::deltaMethod(mvec,"x/y",V,level=.95)
  est90 <- car::deltaMethod(mvec,"x/y",V,level=.90)
  stat <- as.numeric(est[1])/as.numeric(est[2])
  result <- c(Difference = as.numeric(est[1]), SE = as.numeric(est[2]), ll = as.numeric(est[3]), ul = as.numeric(est[4]), ll90 = as.numeric(est90[3]), ul90 = as.numeric(est90[4]),
              pval = 2 * pnorm(-abs(stat)))
  return(result)
}

#acte_fn(df,d2="treatment", dv2="op_inds_1", c_acte_s,w = T)



#### Effects Function ####

effects_fn <- function(the_dat, the_dv,
                       the_lab = c("ATE",
                                   "ACTE:","Choose Information","Avoid Information",
                                   "CATE:",
                                   "Personal Exp","No Personal",
                                   "Major Problem","Not Major",
                                   "College Degree ","No Degree",
                                   "Low Income","High Income",
                                   "White","Non-White",
                                   "Male","Female",
                                   "Democrat","Independent","Republican",
                                   "Liberal","Moderate","Conservative"
                       ),...){
  tmp <- rbind(
    ate_fn(the_dat, dv1=the_dv, c = c("Control","Treatment") ),
    rep(NA,7),
    acte_fn(the_dat, dv2=the_dv, z = c("Control","Selection"),w = T),
    acte_fn(the_dat, dv2=the_dv, z = c("Selection","Treatment"),w = T),
    rep(NA,7),
    ate_fn(the_data = the_dat[the_dat$op_know_personal==1,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$op_know_personal==0,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$op_prob==1,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$op_prob!=1,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$educ>4,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$educ<5,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$income<median(the_dat$income, na.rm=T),], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$income>=median(the_dat$income, na.rm=T),], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$nonwhite01==0,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$nonwhite01==1,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$female01==0,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$female01==1,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$pid<4,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$pid==4,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$pid>4,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$ideology<3,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$ideology==3,], the_dv, c = c("Control","Treatment")),
    ate_fn(the_data = the_dat[the_dat$ideology>3,], the_dv, c = c("Control","Treatment"))
  )
  tmp <- data.frame(tmp)
  #print(dim(tmp))
  tmp$Estimate <- as.character(the_lab)
  tmp$Estimate <- factor(tmp$Estimate,levels=rev(the_lab))
  tmp$Estimand <- factor(c("ATE",
                    rep("ACTE",3),
                    rep("CATE",19)
                    ),
                    levels = c("ATE","ACTE","CATE")
  )
  return(tmp)
}

# effects_fn(df,the_dv = "dv_cause_users")


#### Plot Treatment Effects function ####

plot_effects_fn <- function(d, dv){
  tmp <- effects_fn(the_dat = d, the_dv = dv)
  tmp %>%
    ggplot(aes(Estimate, Difference,col=Estimand))+
    geom_point(aes(shape=Estimand))+
    geom_linerange(aes(ymin=ll,ymax=ul),size=.3)+
    geom_linerange(aes(ymin=ll90,ymax=ul90),size=.6)+
    geom_hline(yintercept = 0,linetype="dashed")+
    theme(axis.text.x = element_text(
      angle = 0, hjust = 1),
      axis.text.y = element_text(face=ifelse(rev(tmp$Estimate)%in%c("ATE","ACTE:","CACTE:","CATE:"),"bold","italic"),
                                 size=ifelse(rev(tmp$Estimate)%in%c("ATE","ACTE:","CACTE:","CATE:"),12,10)
      )
    ) + coord_flip()

}
# plot_effects_fn(df, "dv_cause_providers")

table_effects_fn <- function(d, dv,the_cap="Treatment Effects for Support for Policy Chances"){
  tmp <- effects_fn(d, dv)
  rownames(tmp) <- tmp$Estimate
  tmp <- tmp[-c(2,5),-c(5,6,8,9)]
  tmp <- tmp %>%
    rename(`Pr(<|t|)`="pval")
  tab <- kable(tmp,
               booktabs = TRUE,
               caption = the_cap,
               digits=2,
  align = "l") %>%
  kable_styling(latex_options = c("hold_position",font_size=10)) %>%
    pack_rows("ATE",1,1) %>%
    pack_rows("ACTE",2,3) %>%
    pack_rows("CATE",4,5) %>%
    pack_rows("",6,7) %>%
    pack_rows("",8,9,latex_gap_space = "0.3em") %>%
    pack_rows("",10,11,latex_gap_space = "0.3em") %>%
    pack_rows("",12,13,latex_gap_space = "0.3em") %>%
    pack_rows("",14,15,latex_gap_space = "0.3em") %>%
    pack_rows("",16,18,latex_gap_space = "0.3em") %>%
    pack_rows("",19,21,latex_gap_space = "0.3em")


  return(tab)
  }




# #### Treatment Effects function ####
#
# te_fn <- function(the_dat,the_dv,...){
#   # if(grepl("exp",the_dv)){
#   #   ates <- ates_s
#   #   actes <- actes_s
#   # }else{
#   #   ates <- ates_g
#   #   actes <- actes_g
#   # }
#
#   the_ates <- as.data.frame(do.call(rbind,lapply(ates,function(x)ate_fn(the_dat,dv1=the_dv,c=x,...))))
#   the_actes <- as.data.frame(do.call(rbind,lapply(actes,function(x)acte_fn(the_dat,dv2=the_dv,z=x,w=T,...))))
#   the_s_actes <- s_acte_fn(the_dat,dv3 = the_dv,w2 = F)
#   results<- rbind(the_ates,the_actes,the_s_actes)
#   rownames(results) <- c("ATE","ACTE - Select","ACTE - Avoid", "CACTE - Treatment","CACTE - Alternative")
#   return(results)
#
#
# }
#
# #### Plot Treatment Effects ####
#
# plot_fn <- function(plot_df,to_plot="dv_pca",caption=NULL,...){
#   tmp <- te_fn(the_dat = plot_df,the_dv = to_plot,...)
#   tmp$Effect <- rownames(tmp)
#   tmp$Effect <- factor(tmp$Effect,levels=tmp$Effect)
#   tmp$Estimand <- factor(c("ATE","ACTE","ACTE","CACTE","CACTE"),levels=c("ATE","ACTE","CACTE"))
#   p <- ggplot(tmp, aes(Effect, Estimate,col=Estimand))+
#     geom_point()+
#     geom_errorbar(aes(ymin=ll,ymax=ul),width=.2)+
#     geom_hline(yintercept = 0,linetype="dashed")+
#     theme(axis.text.x = element_text(
#       angle = 90, hjust = 1)
#     )+
#     labs(title=caption)
#   return(p)
#
# }


#### Balace Function ####



balance_fn <- function(the_data, dv1="Y",c=c("Avoid Information","Choose Information"), weights=T,...){
  #frm<-paste(dv1,"~  age+age_sq+is_female+is_nonwhite+is_democrat+is_moderate+is_liberal+open+consc+extra+agree+neur")
  #df[,dv1] <- resid(lm(frm,df,na.action="na.exclude"))
  df <- the_data
  print(dim(df))
  # Weights
  tmp <- as.data.frame(df[df$Treatment%in%c, ])
  #print(dim(tmp))
  if(weights==F){
    mu1 <- with(tmp, mean(tmp[Treatment==c[1], dv1],na.rm=T))
    mu2 <- with(tmp, mean(tmp[Treatment==c[2], dv1],na.rm=T))
  }
  if(weights==T){
    mu1 <- with(tmp, Hmisc::wtd.mean(tmp[Treatment==c[1], dv1],na.rm=T,weights=tmp[Treatment==c[1],"weights"]))
    mu2 <- with(tmp, Hmisc::wtd.mean(tmp[Treatment==c[2], dv1],na.rm=T,weights=tmp[Treatment==c[2],"weights"]))
  }
  diff <- mu1-mu2
  if(weights==F){
    sd1 <- with(tmp, sd(tmp[Treatment==c[1], dv1],na.rm=T))
    sd2 <- with(tmp, sd(tmp[Treatment==c[2], dv1],na.rm=T))
  }
  if(weights==T){
    sd1 <- sqrt(with(tmp, Hmisc::wtd.var(tmp[Treatment==c[1], dv1],na.rm=T,weights = tmp[Treatment==c[1],"weights"])))
    sd2 <- sqrt(with(tmp, Hmisc::wtd.var(tmp[Treatment==c[2], dv1],na.rm=T,weights = tmp[Treatment==c[2],"weights"])))
  }
  n1 <- with(tmp, sum(!is.na(tmp[Treatment==c[1], dv1])*tmp[Treatment==c[1],"weights"]))
  n2 <- with(tmp, sum(!is.na(tmp[Treatment==c[2], dv1])*tmp[Treatment==c[2],"weights"]))
  #print(c(mu1,mu2,n1,n2,sd1,sd2))
  se <- sqrt( sd1^2/n1 + sd2^2/n2)
  the_df <- (sd1^2/n1+sd2^2/n2)^2/((sd1^4)/(n1^2*(n1-1))+ (sd2^4)/(n2^2*(n2-1)))
  ll <- diff - qt(.975,the_df)*se
  ul <- diff + qt(.975,the_df)*se
  ll90 <- diff - qt(.95,the_df)*se
  ul90 <- diff + qt(.95,the_df)*se
  stat <- diff/se
  pval = 2 * pt(-abs(stat),the_df)
  result <- c(Mu1= mu1, Mu2 = mu2, Difference = diff, SE = se, ll = ll, ul = ul,ll90 = ll90, ul90 = ul90, pval = pval, N1=n1,N2=n2)
  return(result)
}


# #### Power Function #####
#
# power_fn <- function(sims,
#                      p_N=1000,
#                      p_sigma = 1,
#                      p_p_treat=.5,
#                      p_prop_select=.5,
#                      p_p_treat_select = c(.25,.5,.25),
#                      p_tau_st = .5,
#                      p_tau_at = -.5,
#                      p_tau_aa = .5,
#                      p_select_effect = 0,
#
#                      ...){
#
#   ate <- rep(NA,sims)
#   acte_s <- rep(NA,sims)
#   acte_a <- rep(NA,sims)
#   acte_at <- rep(NA,sims)
#   acte_aa <- rep(NA,sims)
#
#   sig_ate <- rep(NA,sims)
#   sig_acte_s <- rep(NA,sims)
#   sig_acte_a <- rep(NA,sims)
#   sig_acte_at <- rep(NA,sims)
#   sig_acte_aa <- rep(NA,sims)
#   cor_select <- rep(NA,sims)
#
#
#   for(i in 1:sims){
#     df <- data_fn(N=p_N,sigma= p_sigma ,
#                   p_treat=p_p_treat,
#                   prop_select = p_prop_select,
#                   p_treat_select = p_p_treat_select ,
#                   tau_st= p_tau_st ,
#                   tau_at= p_tau_at ,
#                   tau_aa= p_tau_aa ,
#                   select_effect= p_select_effect )
#     sig_ate[i] <- ate_fn(df,dv1="Y",c=c("Control","Treatment"),weights = T)["pval"]
#     sig_acte_s[i] <- acte_fn(df,z=c_acte_s, w=T)["pval"]
#     sig_acte_a[i] <- acte_fn(df,z=c_acte_a, w=T)["pval"]
#     sig_acte_at[i] <- s_acte_fn(df)[1,"pval"]
#     sig_acte_aa[i] <- s_acte_fn(df)[2,"pval"]
#     cor_select[i] <- cor(df$Y,df$select)
#
#     ate[i] <- ate_fn(df,dv1="Y",c=c("Control","Treatment"))["Difference"]
#     acte_s[i] <- acte_fn(df,z=c_acte_s, w=T)["Difference"]
#     acte_a[i] <- acte_fn(df,z=c_acte_a,w = T)["Difference"]
#     acte_at[i] <- s_acte_fn(df)[1,"Difference"]
#     acte_aa[i] <- s_acte_fn(df)[2,"Difference"]
#   }
#
#   pow_ate <- mean(sig_ate<.05)
#   pow_acte_s <- mean(sig_acte_s<.05)
#   pow_acte_a <- mean(sig_acte_a<.05)
#   pow_acte_at <-  mean(sig_acte_at<.05)
#   pow_acte_aa <- mean(sig_acte_aa<.05)
#   mn_cor_select <- mean(cor_select)
#
#   mn_ate <- mean(ate)
#   mn_acte_s <- mean(acte_s)
#   mn_acte_a <- mean(acte_a)
#   mn_acte_at <-  mean(acte_at)
#   mn_acte_aa <- mean(acte_aa)
#
#
#   res <- rbind(
#     c(mn_ate,
#       mn_acte_s,
#       mn_acte_a,
#       mn_acte_at,
#       mn_acte_aa,mn_cor_select),
#
#     c(pow_ate,
#       pow_acte_s,
#       pow_acte_a,
#       pow_acte_at,
#       pow_acte_aa,
#       NA
#     ))
#
#   return(res)
#
# }
#
# #### Simulate Power ####
#
#
# sim_power_fn <- function(
#                          s_sims = 1000,
#                          s_N=1000,
#                          s_sigma = 1,
#                          s_p_treat=.5,
#                          s_prop_select=.5,
#                          s_p_treat_select = c(.25,.5,.25),
#                          s_tau_st = .5,
#                          s_tau_at = -.5,
#                          s_tau_aa = .5,
#                          s_select_effect = 0
#
#
#
# ) {
#   # Create matrix to store values
#   power_mat <- matrix(NA,nrow=5, ncol = length(s_tau_st),
#                        dimnames = list(c("ATE","ACTE-Select","ACTE-Avoid","CACTE-Treatment","CACTE-Alternative"),
#                                        s_tau_st
#                        )
#                       )
#   bias_mat <- matrix(NA,nrow=5, ncol = length(s_tau_st),
#                       dimnames = list(c("ATE","ACTE-Select","ACTE-Avoid","CACTE-Treatment","CACTE-Alternative"),
#                                       s_tau_st
#                       )
#   )
#   tmp <- c()
#   ave_cor <- c()
#   tmp_df <- data.frame(Estimate=NULL,Tau=NULL)
#   df <- data.frame(Estimate= NULL,
#                    Tau_Select = NULL,
#                    Tau_Avoid = NULL ,
#                    Tau_Alt = NULL ,
#                    Power = NULL )
#   # Loop over possible values
#   for(i in 1:length(s_tau_st)){
#     tmp <- power_fn(sims = s_sims,
#                                p_N = s_N,
#                                p_sigma  = s_sigma,
#                                p_p_treat = s_p_treat,
#                                p_prop_select = s_prop_select,
#                                p_p_treat_select  = s_p_treat_select,
#                                p_tau_st  = s_tau_st[i],
#                                p_tau_at  = s_tau_at[i],
#                                p_tau_aa  = s_tau_aa[i],
#                                p_select_effect  = s_select_effect
#     )
#     ave_cor[i] <- tmp[1,6]
#     power_mat[,i] <- tmp[2,1:5]
#     tmp2 <- data.frame(Estimate=c("ATE","ACTE-Select","ACTE-Avoid","CACTE-Treatment","CACTE-Alternative"),
#                        Type = c("ATE","ACTE","ACTE","CACTE","CACTE"),
#                   Tau_Select = rep(s_tau_st[i],5),
#                   Tau_Avoid = rep(s_tau_at[i],5),
#                   Tau_Alt = rep(s_tau_aa[i],5),
#                   Power = tmp[2,1:5]
#                   )
#     df <- rbind(tmp2,df)
#   }
#   return(list(df,power_mat,ave_cor))
#
# }
#
# power_fn(sims = 4)
#
# #### Plot Power Simulations ####
#
#
# plot_power_sim_fn <- function(
#   p_s_sims = 1000,
#   p_s_N=1000,
#   p_s_sigma = 1,
#   p_s_p_treat=.4,
#   p_s_prop_select=.5,
#   p_s_p_treat_select = c(.25,.5,.25),
#   p_s_tau_st = .5,
#   p_s_tau_at = -.5,
#   p_s_tau_aa = .5,
#   p_s_select_effect = 0,
#   lab_effects = "Effects = Equal & offsetting"
#
# ){
#   pow <- sim_power_fn(
#     s_sims = p_s_sims,
#     s_N= p_s_N,
#     s_sigma = p_s_sigma,
#     s_p_treat= p_s_p_treat,
#     s_prop_select= p_s_prop_select,
#     s_p_treat_select = p_s_p_treat_select,
#     s_tau_st = p_s_tau_st,
#     s_tau_at = p_s_tau_at,
#     s_tau_aa = p_s_tau_aa,
#     s_select_effect = p_s_select_effect
#
#
#   )
#   rhos.min <- round(range(pow[[3]]),2) [1]
#   rhos.max <- round(range(pow[[3]]),2) [2]
#
#   p <- pow[[1]] %>%
#     ggplot(aes(Tau_Select,Power, col=Estimate,linetype=Type))+
#     geom_line()+
#     ylim(0,1.05) +
#     xlim(0,.7) +
#     geom_hline(yintercept = .8,linetype = "dashed",col="grey")+
#     xlab(expression(tau[Select]))+
#     ylab("Power (Probability of Statistical Signficance)")+
#     annotate(geom = "text",
#              hjust = 0,
#              y = 1.05,
#              x = 0,
#              label = paste("Simulations =",scales::comma(p_s_sims))
#     )+
#     annotate(geom = "text",
#              hjust = 0,
#              y = 1,
#              x = 0,
#              label = lab_effects
#     )+
#     annotate(geom = "text",
#              hjust = 0,
#              y = .95,
#              x = 0,
#              label =paste("Prop Select =", round(p_s_prop_select,2))
#     ) +
#     annotate(geom = "text",
#              hjust = 0,
#              y = .90,
#              x = 0,
#              label = paste("Cor = [",rhos.min,", ",rhos.max,"]",sep="")
#     )
#   tab <- kable(pow[[2]],
#                caption = "Power Analysis",
#                format = "latex",
#                booktabs=T,
#                linesep = "",
#                digits=2) %>%
#     add_header_above(c("","Hypothesized Effect Among Selectors"=dim(pow[[2]])[2])) %>%
#     kable_styling(latex_options = c("hold_position",font_size=10))
#
#   return(list(p,tab))
#
#
# }
#
#
# ### Sim Bias ###
#
# bias_fn <- function(sims,
#                      p_N=1000,
#                      p_sigma = 1,
#                      p_p_treat=.5,
#                      p_prop_select=.5,
#                      p_p_treat_select = c(.25,.5,.25),
#                      p_tau_st = .5,
#                      p_tau_at = -.5,
#                      p_tau_aa = .5,
#                      p_select_effect = 0,
#
#                      ...){
#
#   ate <- rep(NA,sims)
#   acte_s <- rep(NA,sims)
#   acte_a <- rep(NA,sims)
#   acte_at <- rep(NA,sims)
#   acte_aa <- rep(NA,sims)
#
#   sig_ate <- rep(NA,sims)
#   sig_acte_s <- rep(NA,sims)
#   sig_acte_a <- rep(NA,sims)
#   sig_acte_at <- rep(NA,sims)
#   sig_acte_aa <- rep(NA,sims)
#   cor_select <- rep(NA,sims)
#
#
#   for(i in 1:sims){
#     df <- data_fn(N=p_N,sigma= p_sigma ,
#                   p_treat=p_p_treat,
#                   prop_select = p_prop_select,
#                   p_treat_select = p_p_treat_select ,
#                   tau_st= p_tau_st ,
#                   tau_at= p_tau_at ,
#                   tau_aa= p_tau_aa ,
#                   select_effect= p_select_effect )
#     sig_ate[i] <- ate_fn(df,dv1="Y",c=c("Control","Treatment"),weights = T)["pval"]
#     sig_acte_s[i] <- acte_fn(df,z=c_acte_s, w=T)["pval"]
#     sig_acte_a[i] <- acte_fn(df,z=c_acte_a, w=T)["pval"]
#     sig_acte_at[i] <- s_acte_fn(df)[1,"pval"]
#     sig_acte_aa[i] <- s_acte_fn(df)[2,"pval"]
#     cor_select[i] <- cor(df$Y,df$select)
#
#     ate[i] <- ate_fn(df,dv1="Y",c=c("Control","Treatment"))["Difference"]
#     acte_s[i] <- acte_fn(df,z=c_acte_s, w=T)["Difference"]
#     acte_a[i] <- acte_fn(df,z=c_acte_a,w = T)["Difference"]
#     acte_at[i] <- s_acte_fn(df)[1,"Difference"]
#     acte_aa[i] <- s_acte_fn(df)[2,"Difference"]
#   }
#
#   pow_ate <- mean(sig_ate<.05)
#   pow_acte_s <- mean(sig_acte_s<.05)
#   pow_acte_a <- mean(sig_acte_a<.05)
#   pow_acte_at <-  mean(sig_acte_at<.05)
#   pow_acte_aa <- mean(sig_acte_aa<.05)
#   mn_cor_select <- mean(cor_select)
#
#   mn_ate <- mean(ate)
#   mn_acte_s <- mean(acte_s)
#   mn_acte_a <- mean(acte_a)
#   mn_acte_at <-  mean(acte_at)
#   mn_acte_aa <- mean(acte_aa)
#
#
#   res <- cbind(
#     ate,
#     acte_s,
#     acte_a,
#     acte_at,
#     acte_aa)
#   return(res)
#
# }

# tmp <- bias_fn(500,p_prop_select = .8,p_tau_st = .5,p_tau_at = -.3,p_select_effect = -2)
#
# plot(density(tmp[,2]))
# abline(v=0.5)
# abline(v=mean(tmp[,2]),col="red")
#
# plot(density(tmp[,3]))
# abline(v=-.3)
# abline(v=mean(tmp[,3]),col="red")
#
# mean(tmp[,4]<0)
#
# plot(density(tmp[,4]))
# abline(v=-3)
# abline(v=mean(tmp[,4]),col="red")
#
# plot(density(tmp[,5]))
# abline(v=.5)
# abline(v=mean(tmp[,5]),col="red")
#
#
# # plot_power_sim_fn(p_s_sims = 2,
# #                   p_s_tau_st = seq(.1,.7,by=.1),
# #                   p_s_tau_at = seq(-.1,-.7,by=-.1),
# #                   p_s_tau_aa = seq(.1,.7,by=.1)
# #
# #                   )
